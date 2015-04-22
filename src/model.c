#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "model.h"
#include "vector.h"
#include "sterics.h"

/**
 * Allocate memory for a model structure.
 *
 * \return A newly-allocated model structure or NULL if out of memory.
 */
struct model *model_alloc(){
    struct model *m = malloc(sizeof(struct model));
    if(!m)
        return NULL;

    m->num_linear_springs = 0;
    m->num_torsion_springs = 0;
    m->num_residues = 0;
    m->residues = NULL;
    m->linear_springs = NULL;
    m->torsion_springs = NULL;
    m->time = 0;
    m->until = 0;
    m->timestep = 0.1;
    m->synth_time = 100;
    m->drag_coefficient = -0.1;
    m->shield_drag = false;
    m->steric_grid = NULL;
    m->use_sterics = false;
    m->use_water = false;
    m->max_synth_angle = DEFAULT_MAX_SYNTH_ANGLE;
    m->fix = false;
    m->threestate = false;
    m->do_synthesis = true;
    return m;
}

/**
 * Free memory for a model structure, including the residues, atoms and springs.
 */
void model_free(struct model *m){
    for(size_t i=0; i < m->num_residues; i++)
        free(m->residues[i].atoms);
    free(m->residues);
    free(m->linear_springs);
    free(m->torsion_springs);
    free(m);
}

void model_accumulate_forces(struct model *m){
    struct vector force;
    struct vector tmp;

    struct residue *residues = m->residues;
    struct linear_spring *linear_springs = m->linear_springs;
    struct torsion_spring *torsion_springs = m->torsion_springs;

    //Begin by zeroing out any existing forces
    #pragma omp parallel for shared(residues)
    for(size_t i=0; i < m->num_residues; i++){
        for(size_t j=0; j < residues[i].num_atoms; j++)
            vector_zero(&residues[i].atoms[j].force);
    }

    //Then go through all springs and accumulate forces on the residues
    #pragma omp parallel for shared(linear_springs)
    for(size_t i=0; i < m->num_linear_springs; i++){
        struct linear_spring *s = &linear_springs[i];

        if(s->a->synthesised && s->b->synthesised){
            if(!s->a->fixed){
                linear_spring_force(&force, s, A);
                vadd_to(&s->a->force, &force);
            }

            if(!s->b->fixed){
                linear_spring_force(&force, s, B);
                vadd_to(&s->b->force, &force);
            }
        }
    }

    //Torsion springs
    #pragma omp parallel for shared(torsion_springs)
    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring *s = &torsion_springs[i];

        if(s->a1->synthesised
                && s->a2->synthesised
                && s->a3->synthesised
                && s->a4->synthesised){

            if(!s->a1->fixed){
                torsion_spring_force(&force, s, R1);
                vadd_to(&s->a1->force, &force);
            }

            if(!s->a4->fixed){
                torsion_spring_force(&force, s, R4);
                vadd_to(&s->a4->force, &force);
            }
        }
    }

    //If we're not using the fancy drag force, apply the drag force now.
    if(!m->shield_drag){
        for(size_t i=0; i < m->num_residues; i++){
            for(size_t j=0; j < residues[i].num_atoms; j++){
                vector_copy_to(&tmp, &residues[i].atoms[j].velocity);
                vmul_by(&tmp, m->drag_coefficient);
                vadd_to(&residues[i].atoms[j].force, &tmp);
            }
        }
    }

    //Steric, water and drag forces
    if(m->steric_grid){
        steric_grid_update(m->steric_grid, m);

        if(m->use_sterics)
            steric_grid_forces(m->steric_grid, m);
        if(m->use_water)
            water_force(m, m->steric_grid);
        if(m->shield_drag)
            drag_force(m, m->steric_grid);
    }
}


const char *atom_fmt   = "ATOM  %5d  %-3s %-3s  %4d%1s   %8.3f%8.3f%8.3f\n";
const char *conect_fmt = "CONECT% 5d% 5d\n";

int model_pdb(FILE *out, const struct model *m, bool conect){

    int bytes_written = 0;
    for(size_t i=0; i < m->num_residues; i++){
        struct residue r = m->residues[i];

        for(size_t j=0; j < r.num_atoms; j++){
            struct atom a = r.atoms[j];
            if(a.synthesised){
                int res = fprintf(out, atom_fmt, a.id, a.name,
                        r.name,
                        r.id,
                        " ",
                        a.position.c[0],
                        a.position.c[1],
                        a.position.c[2]);

                if(res < 0)
                    return res;
                bytes_written += res;
            }
        }
    }
    if(conect){
        for(size_t i=0; i < m->num_linear_springs; i++){
            struct linear_spring s = m->linear_springs[i];
            if(linear_spring_active(&s) && s.a->synthesised && s.b->synthesised){
                int res = fprintf(out, conect_fmt, s.a->id, s.b->id);

                if(res < 0)
                    return res;
                bytes_written += res;
            }
        }
    }
    return bytes_written;
}

/**
 * Generate a model at time src->time.
 *
 * \param[out] dst State will be placed in here.
 * \param[in]  src Full model.
 *
 * \details This function copies the state of model src to dst such that only
 * residues that have been synthesised in the past are included. This is done
 * by changing the num_residues parameter of src and using the same array of
 * residues, so don't free any members of dst!
 */
void model_synth(struct model *dst, const struct model *src){
    memcpy(dst, src, sizeof(*src));

    if(!src->do_synthesis)
        dst->num_residues = src->num_residues;
    else
        dst->num_residues = src->time / src->synth_time;

    if(dst->num_residues > src->num_residues)
        dst->num_residues = src->num_residues;

    for(size_t i=0; i < dst->num_residues; i++){
        struct residue *prev  = (i >= 1) ? &dst->residues[i-1] : NULL;
        struct residue *prev2 = (i >= 2) ? &dst->residues[i-2] : NULL;
        if(!dst->residues[i].synthesised){
            residue_synth(&dst->residues[i], prev, prev2, src->max_synth_angle);
            if(dst->fix && prev){
                for(size_t j=0; j < prev->num_atoms; j++){
                    prev->atoms[j].fixed = true;
                    vector_zero(&prev->atoms[j].velocity);
                }
            }
        }

    }
}

/**
 * Add a residue to the model.
 *
 * This reallocates the residues array and increases num_residues.
 *
 * \param id ID of the residue to add.
 */
struct residue * model_push_residue(struct model *m, int id){
    m->residues = realloc(m->residues,
            sizeof(struct residue) * (m->num_residues + 1));
    residue_init(&m->residues[m->num_residues], id);
    m->num_residues++;
    return &m->residues[m->num_residues - 1];
}
