#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "model.h"
#include "vector.h"
#include "sterics.h"
#include "rama.h"
#include "linear_spring.h"
#include "bond_angle.h"
#include "rama.h"
#include "torsion_spring.h"

static void add_torsion_force(struct torsion_spring *spring);

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
    m->num_rama_constraints = 0;
    m->num_bond_angles = 0;
    m->num_residues = 0;
    m->num_atoms = 0;
    m->residues = NULL;
    m->atoms = NULL;
    m->linear_springs = NULL;
    m->torsion_springs = NULL;
    m->bond_angles = NULL;
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
    m->debug = NULL;
    return m;
}

/**
 * Free memory for a model structure, including the residues, atoms and springs.
 */
void model_free(struct model *m){
    free(m->atoms);
    free(m->residues);
    free(m->linear_springs);
    free(m->torsion_springs);
    free(m->bond_angles);
    free(m->rama_constraints);
    free(m);
}

void model_accumulate_forces(struct model *m){
    struct vector force;
    struct vector tmp;

    struct linear_spring *linear_springs = m->linear_springs;
    struct torsion_spring *torsion_springs = m->torsion_springs;
    struct bond_angle_spring *bond_angles = m->bond_angles;

    //Begin by zeroing out any existing forces
    for(size_t i=0; i < m->num_atoms; i++)
        vector_zero(&m->atoms[i].force);

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

            //Print debug information
            if(m->debug && m->debug->linear
                    && (int)(m->time / m->debug->interval) > m->debug->nprinted){
                struct vector force_a, force_b;
                struct vector displacement;
                linear_spring_force(&force_a, s, A);
                linear_spring_force(&force_b, s, B);
                vsub(&displacement, &s->a->position, &s->b->position);

                fprintf(m->debug->linear, DEBUG_LINEAR_FMT,
                        m->time,
                        s->a->id, s->a->name,
                        s->b->id, s->b->name,
                        (s->enabled ? "enabled" : "disabled"),
                        s->distance, vmag(&displacement), 
                        force_a.c[0], force_a.c[1], force_a.c[2],
                        force_b.c[0], force_b.c[1], force_b.c[2]);
            }
        }
    }

    //Torsion springs
    #pragma omp parallel for shared(torsion_springs)
    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring *s = &torsion_springs[i];
        add_torsion_force(s);

        if(m->debug && m->debug->torsion
            && (int)(m->time / m->debug->interval) > m->debug->nprinted){
            if(!s->a1->synthesised || !s->a2->synthesised
                    || !s->a3->synthesised || !s->a4->synthesised)
                continue;



            struct vector spring_forces[4];
            torsion_spring_force_new(
                    &spring_forces[0],
                    &spring_forces[1],
                    &spring_forces[2],
                    &spring_forces[3],
                    s);
                fprintf(m->debug->torsion, DEBUG_TORSION_FMT,
                        m->time,
                        s->a1->id, s->a1->name,
                        s->a2->id, s->a2->name,
                        s->a3->id, s->a3->name,
                        s->a4->id, s->a4->name,
                        (s->enabled ? "enabled" : "disabled"),
                        s->angle, torsion_spring_angle(s),
                        spring_forces[0].c[0],
                        spring_forces[0].c[1],
                        spring_forces[0].c[2],
                        spring_forces[1].c[0],
                        spring_forces[1].c[1],
                        spring_forces[1].c[2],
                        spring_forces[2].c[0],
                        spring_forces[2].c[1],
                        spring_forces[2].c[2],
                        spring_forces[3].c[0],
                        spring_forces[3].c[1],
                        spring_forces[3].c[2]);
        }
    }

    //Ramachandran constraints
    for(size_t i=0; i < m->num_rama_constraints; i++){
        struct rama_constraint *rama = &m->rama_constraints[i];
        if(rama_is_synthesised(rama)){
            rama_get_closest(rama);
            if(rama->enabled){
                add_torsion_force(rama->phi);
                add_torsion_force(rama->psi);
            }
        }
    }

    //Bond angle constraints
    #pragma omp parallel for shared(torsion_springs)
    for(size_t i=0; i < m->num_bond_angles; i++){
        struct bond_angle_spring *s = &bond_angles[i];

        if(s->a1->synthesised
                && s->a2->synthesised
                && s->a3->synthesised){

            struct vector spring_forces[3];
            bond_angle_force(
                    &spring_forces[0],
                    &spring_forces[1],
                    &spring_forces[2],
                    s);

            if(!s->a1->fixed)
                vadd_to(&s->a1->force, &spring_forces[0]);
            if(!s->a2->fixed)
                vadd_to(&s->a2->force, &spring_forces[1]);
            if(!s->a3->fixed)
                vadd_to(&s->a3->force, &spring_forces[2]);

            //Print debug information
            if(m->debug && m->debug->angle
                && (int)(m->time / m->debug->interval) > m->debug->nprinted){
                fprintf(m->debug->angle, DEBUG_ANGLE_FMT,
                        m->time,
                        s->a1->id, s->a1->name,
                        s->a2->id, s->a2->name,
                        s->a3->id, s->a3->name,
                        (s->enabled ? "enabled" : "disabled"),
                        s->angle, bond_angle_angle(s),
                        spring_forces[0].c[0],
                        spring_forces[0].c[1],
                        spring_forces[0].c[2],
                        spring_forces[1].c[0],
                        spring_forces[1].c[1],
                        spring_forces[1].c[2],
                        spring_forces[2].c[0],
                        spring_forces[2].c[1],
                        spring_forces[2].c[2]);
            }
        }
    }

    //If we're not using the fancy drag force, apply the drag force now.
    if(!m->shield_drag){
        for(size_t i=0; i < m->num_atoms; i++){
            vector_copy_to(&tmp, &m->atoms[i].velocity);
            vmul_by(&tmp, m->drag_coefficient);
            vadd_to(&m->atoms[i].force, &tmp);
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

    if(m->debug && m->debug->angle
        && (int)(m->time / m->debug->interval) > m->debug->nprinted){
        m->debug->nprinted++;
    }
}

static void add_torsion_force(struct torsion_spring *s){
    if(s->a1->synthesised
            && s->a2->synthesised
            && s->a3->synthesised
            && s->a4->synthesised){

        struct vector spring_forces[4];
        torsion_spring_force_new(
                &spring_forces[0],
                &spring_forces[1],
                &spring_forces[2],
                &spring_forces[3],
                s);

        if(!s->a1->fixed)
            vadd_to(&s->a1->force, &spring_forces[0]);
        if(!s->a2->fixed)
            vadd_to(&s->a2->force, &spring_forces[1]);
        if(!s->a3->fixed)
            vadd_to(&s->a3->force, &spring_forces[2]);
        if(!s->a4->fixed)
            vadd_to(&s->a4->force, &spring_forces[3]);
    }
}



const char *atom_fmt   = "ATOM  %5d  %-3s %-3s  %4d%1s   %8.3f%8.3f%8.3f\n";
const char *conect_fmt = "CONECT% 5d% 5d\n";

int model_pdb(FILE *out, const struct model *m, bool conect){

    int bytes_written = 0;
    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a    = &m->atoms[i];
        struct residue *r = &m->residues[a->residue_idx];

        if(a->synthesised){
            fprintf(out, atom_fmt, a->id, a->name,
                    r->name,
                    r->id,
                    " ",
                    a->position.c[0],
                    a->position.c[1],
                    a->position.c[2]);

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


/** Synthesise (i.e. set position and synthesised flag) a new atom.
 *
 * Here, an atom is "synthesised" by placing all atoms into the system and
 * setting the "synthesised" property of the residue and all atoms to true.
 * The way that atoms are synthesised depends on both the number of atoms
 * already in the system and the type of atom being synthesised.
 *
 * - If this atom is the first atom being synthesised, just place it at the
 *   origin.
 *
 * - If the second atom is being placed and it is a backbone atom, move it
 *   along the z-axis, plus a random azimuthal angle, according to the steric
 *   radii of the two atoms (i.e. place them so that they are not quite
 *   touching.) If it is not a backbone atom, pick a direction in the x-y plane
 *   for it using the same distance rules.
 *
 * - For any other atoms, define the backbone axis as the direction taken by
 *   the previous two backbone atoms. Backbone atoms get placed along this
 *   vector, and non-backbone atoms get placed around this vector.
 */
void model_synth_atom(const struct model *m, size_t idx, double max_angle){
    //We are going to be synthesising this atom:
    struct atom *a = &m->atoms[idx];

    //Try and get the previous two backbone atoms
    struct atom *prev1 = NULL;
    struct atom *prev2 = NULL;

    for(size_t i=idx - 1; idx > 0 && i >= 0 && !prev1 && !prev2; i--){
        if(m->atoms[i].backbone){
            if(prev1)
                prev2 = &m->atoms[i];
            else
                prev1 = &m->atoms[i];
        }
    }

    a->synthesised = true;
    if(!prev1 && !prev2){
        //If this is the first atom, just plonk it down
        vector_zero(&a->position);
    }else if(prev1 && !prev2){
        //If this is the second atom, just plonk it down near the z-axis for a
        //backbone atom, or on the x-y plane for a non-backbone atom.

        //First, find the required distance between this and the previous atom
        double separation = a->radius + prev1->radius;

        struct vector unit_offset;
        if(a->backbone){
            //Get a unit vector near the z-axis
            vector_rand(&unit_offset, 0, max_angle / 180 * M_PI);
        }else{
            //Get a unit vector along the x-y axis
            vector_rand(&unit_offset,
                    (90 - max_angle) / 180 * M_PI,
                    (90 + max_angle) / 180 * M_PI);
        }

        //Multiply by the distance
        vmul_by(&unit_offset, separation);

        //Add to the previous atom's coordinates
        vadd_to(&unit_offset, &prev1->position);

        //Copy into the new atom's coords
        vector_copy_to(&a->position, &unit_offset);
    }else{

        //We want to do the same as before, but then we want to rotate it such
        //that the vector between the previous two atoms is the new z-axis.
        double separation = a->radius + prev1->radius;
        struct vector unit_offset;
        if(a->backbone){
            vector_rand(&unit_offset, 0, max_angle / 180 * M_PI);
        }else{
            vector_rand(&unit_offset,
                    (90 - max_angle) / 180 * M_PI,
                    (90 + max_angle) / 180 * M_PI);
        }
        vmul_by(&unit_offset, separation);
        vadd_to(&unit_offset, &prev1->position);

        //So we have our vector, but we need to rotate our coordinate system.
        //To do this we rotate it with the same angle and axis that the
        //prev1-prev2 displacement vector has from the z-axis.
        struct vector z = {{0, 0, 1}};
        struct vector displacement;
        struct vector rot_axis;
        double angle;
        vsub(&displacement, &prev1->position, &prev2->position);
        vcross(&rot_axis, &displacement, &z);
        angle = acos(vdot(&displacement, &z) / vmag(&displacement));

        vrot_axis(&unit_offset, &rot_axis, &unit_offset, -angle);
        vadd(&a->position, &unit_offset, &prev1->position);
    }
}

