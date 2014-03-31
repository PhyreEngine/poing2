#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "model.h"
#include "vector.h"

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
    m->timestep = 0.1;
    m->synth_time = 100;
    m->drag_coefficient = 0.1;
    return m;
}

/**
 * Free memory for a model structure, including the residues and springs.
 */
void model_free(struct model *m){
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
    for(size_t i=0; i < m->num_residues; i++)
        vector_zero(&residues[i].force);

    //Then go through all springs and accumulate forces on the residues
    #pragma omp parallel for shared(linear_springs)
    for(size_t i=0; i < m->num_linear_springs; i++){
        struct linear_spring s = linear_springs[i];

        if(s.a->synthesised && s.b->synthesised){
            linear_spring_force(&force, &s, A);
            vadd_to(&s.a->force, &force);

            linear_spring_force(&force, &s, B);
            vadd_to(&s.b->force, &force);
        }
    }

    #pragma omp parallel for shared(torsion_springs)
    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring s = torsion_springs[i];

        if(s.r1->synthesised
                && s.r2->synthesised
                && s.r3->synthesised
                && s.r4->synthesised){

            torsion_spring_force(&force, &s, R1);
            vadd_to(&s.r1->force, &force);

            torsion_spring_force(&force, &s, R4);
            vadd_to(&s.r4->force, &force);
        }
    }

    #pragma omp parallel for shared(residues)
    for(size_t i=0; i < m->num_residues; i++){
        vector_copy_to(&tmp, &residues[i].velocity);
        vmul_by(&tmp, -1);
        vadd_to(&residues[i].force, &tmp);
    }
}


const char *atom_fmt   = "ATOM  %5d %4s %-3s  %4d%1s   %8.3f%8.3f%8.3f\n";
const char *conect_fmt = "CONECT% 5d% 5d\n";

char * model_pdb(const struct model *m, bool conect){
    char *buffer = malloc((m->num_residues + m->num_linear_springs) * 80);
    buffer[0] = '\0';

    char line[80];
    line[0] = '\0';
    for(size_t i=0; i < m->num_residues; i++){
        struct residue r = m->residues[i];
        if(r.synthesised){
            sprintf(line, atom_fmt, r.id, " CA ",
                    r.aa->threeletter,
                    r.id,
                    " ",
                    r.position.c[0],
                    r.position.c[1],
                    r.position.c[2]);
        }
        strcat(buffer, line);
    }
    if(conect){
        for(size_t i=0; i < m->num_linear_springs; i++){
            struct linear_spring s = m->linear_springs[i];
            if(s.a->synthesised && s.b->synthesised){
                sprintf(line, conect_fmt, s.a->id, s.b->id);
                strcat(buffer, line);
            }
        }
    }
    return buffer;
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
    dst->num_residues = src->time / src->synth_time;

    if(dst->num_residues > src->num_residues)
        dst->num_residues = src->num_residues;

    for(size_t i=0; i < dst->num_residues; i++){

        /*
         * Find the position to plop down any newly synthesised residues.
         *
         * For the first residue, just leave it where it is.
         *
         * For the second residue, plop it down in some arbitrary direction (in
         * this case, we just use the z direction).
         *
         * For other residues, estimate the position by just adding the vector
         * connecting the previous two points to the position of the previous
         * point.
         *
         */

        if(!dst->residues[i].synthesised){
            if(i == 1){
                //TODO: Set up some kind of initial bond length
                struct vector tmp;
                vector_fill(&tmp, 0, 0, 1);
                for(size_t i = 0; i < N; i++)
                    tmp.c[i] += ((double)rand()) / RAND_MAX * 0.1;

                vadd(&dst->residues[i].position,
                        &tmp,
                        &dst->residues[i-1].position);
            }else if(i >= 2){
                struct vector randv;
                for(size_t i = 0; i < N; i++)
                    randv.c[i] = ((double)rand()) / RAND_MAX * 0.1;
                vector_copy_to(&dst->residues[i].position,
                        &dst->residues[i-1].position);
                vadd_to(&dst->residues[i].position,
                        &dst->residues[i-1].position);
                vsub_to(&dst->residues[i].position,
                        &dst->residues[i-2].position);
                vsub_to(&dst->residues[i].position,
                        &randv);
            }
        }
        dst->residues[i].synthesised = true;
    }
}
