#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "model.h"
#include "vector.h"
#include "sterics.h"
#include "rama.h"
#include "linear_spring.h"
#include "bond_angle.h"
#include "rama.h"
#include "torsion_spring.h"
#include "debug.h"

#ifdef HAVE_CLOCK_GETTIME
#include "profile.h"
#endif

static void add_torsion_force(struct torsion_spring *spring);
static double model_get_separation(
        const struct model *restrict m,
        const struct atom *restrict a,
        const struct atom *restrict b,
        struct atom **restrict place_near);

static void apply_spring_force(struct model *m);
static void apply_torsion_force(struct model *m);
static void apply_rama_force(struct model *m);
static void apply_angle_force(struct model *m);
static void apply_drag_force(struct model *m);
static void profile(struct model *m, const char *msg);

static void model_move_along_vector(struct model *m, double alpha,
        struct vector *r, struct vector *p);

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
    m->num_constraints = 0;
    m->num_residues = 0;
    m->num_atoms = 0;
    m->residues = NULL;
    m->atoms = NULL;
    m->linear_springs = NULL;
    m->torsion_springs = NULL;
    m->bond_angles = NULL;
    m->constraints = NULL;
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
    m->fix_before = -1;
    m->record_time = m->timestep * 10;
    m->max_jitter = 0.01;
    m->profiler = NULL;
    m->bond_map = NULL;
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
    free(m->constraints);
    free(m);
}

void model_accumulate_forces(struct model *m){
    //Begin by zeroing out any existing forces
    for(size_t i=0; i < m->num_atoms; i++)
        vector_zero(&m->atoms[i].force);

    #ifdef HAVE_CLOCK_GETTIME
    if(m->profiler)
        profile_start(m->profiler);
    #endif

    apply_spring_force(m);
    profile(m, "linear");

    apply_torsion_force(m);
    profile(m, "torsion");

    apply_rama_force(m);
    profile(m, "rama");

    apply_angle_force(m);
    profile(m, "angle");

    apply_drag_force(m);
    profile(m, "drag");

    //Steric, water and drag forces
    if(m->steric_grid){
        steric_grid_update(m->steric_grid, m);
        steric_grid_build_ilists(m->steric_grid, m);
        profile(m, "steric grid update");

        if(m->use_sterics){
            steric_grid_forces(m->steric_grid, m);
            profile(m, "steric force");
        }if(m->use_water){
            water_force(m, m->steric_grid);
            profile(m, "water force");
        }if(m->shield_drag){
            drag_force(m, m->steric_grid);
            profile(m, "shielded drag");
        }
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

int model_pdb(FILE *out, const struct model *m, bool conect, int *n){
    int bytes_written = 0;

    bytes_written += fprintf(out, "MODEL     %lu\n", ++(*n));
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
        for(size_t i=0; i < m->num_constraints; i++){
            struct constraint *s = &m->constraints[i];
            struct atom *a = &m->atoms[s->a];
            struct atom *b = &m->atoms[s->b];
            if(a->synthesised && b->synthesised){
                int res = fprintf(out, conect_fmt, a->id, b->id);
                if(res < 0)
                    return res;
                bytes_written += res;
            }
        }
    }
    bytes_written += fprintf(out, "ENDMDL\n");
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

    for(int i=idx - 1; idx > 0 && i >= 0 && !prev2; i--){
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
        struct atom *place_near = prev1;
        double separation = model_get_separation(m, a, prev1, &place_near);

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
        vadd_to(&unit_offset, &place_near->position);

        //Copy into the new atom's coords
        vector_copy_to(&a->position, &unit_offset);
    }else{

        //We want to do the same as before, but then we want to rotate it such
        //that the vector between the previous two atoms is the new z-axis.
        struct atom *place_near = prev1;
        double separation = model_get_separation(m, a, prev1, &place_near);
        struct vector unit_offset;
        if(a->backbone){
            vector_rand(&unit_offset, 0, max_angle / 180 * M_PI);
        }else{
            vector_rand(&unit_offset,
                    (90 - max_angle) / 180 * M_PI,
                    (90 + max_angle) / 180 * M_PI);
        }
        vmul_by(&unit_offset, separation);

        //So we have our vector, but we need to rotate our coordinate system.
        //To do this we rotate it with the same angle and axis that the
        //prev1-prev2 displacement vector has from the z-axis.
        struct vector z = {{0, 0, 1}};
        struct vector displacement;
        struct vector rot_axis;
        double angle;
        vsub(&displacement, &prev1->position, &prev2->position);
        vcross(&rot_axis, &displacement, &z);
        vdiv_by(&rot_axis, vmag(&rot_axis));
        angle = acos(vdot(&displacement, &z) / vmag(&displacement));

        struct vector vout;
        vrot_axis(&vout, &rot_axis, &unit_offset, -angle);
        vadd(&a->position, &vout, &place_near->position);
    }
}

//If the model contains a constraint with these atoms, use that distance as
//the desired distance. We just do a linear search because this function
//will be called relatively infrequently.
//
//If no constraint is found, return the sum of the atom radii.
double model_get_separation(
        const struct model *restrict m,
        const struct atom *restrict a,
        const struct atom *restrict b,
        struct atom **restrict place_near){

    for(size_t i=0; i < m->num_constraints; i++){
        struct atom *c_a = &m->atoms[m->constraints[i].a];
        struct atom *c_b = &m->atoms[m->constraints[i].b];
        if(a == c_a && c_b->synthesised){
            *place_near = c_b;
            return m->constraints[i].distance;
        }else if(a == c_b && c_a->synthesised){
            *place_near = c_a;
            return m->constraints[i].distance;
        }
    }
    return a->radius + b->radius;
}

void apply_spring_force(struct model *m){
    struct vector force1, force2;
    struct linear_spring *linear_springs = m->linear_springs;

    //Then go through all springs and accumulate forces on the residues
    #ifdef HAVE_OPENMP
    #pragma omp parallel for shared(linear_springs)
    #endif
    for(size_t i=0; i < m->num_linear_springs; i++){
        struct linear_spring *s = &linear_springs[i];

        if(s->a->synthesised && s->b->synthesised){
            if(!s->a->fixed || !s->b->fixed)
                linear_spring_force(&force1, &force2, s);

            if(!s->a->fixed)
                vadd_to(&s->a->force, &force1);

            if(!s->b->fixed)
                vadd_to(&s->b->force, &force2);

            //Print debug information
            if(m->debug)
                debug_linear(m, s);
        }
    }
}

void apply_torsion_force(struct model *m){
    struct torsion_spring *torsion_springs = m->torsion_springs;

    //Torsion springs
    #ifdef HAVE_OPENMP
    #pragma omp parallel for shared(torsion_springs)
    #endif
    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring *s = &torsion_springs[i];

        if(s->a1->fixed && s->a2->fixed && s->a3->fixed && s->a4->fixed)
            continue;

        add_torsion_force(s);

        if(m->debug)
            debug_torsion(m, s);
    }
}

void apply_rama_force(struct model *m){
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

}

void apply_angle_force(struct model *m){
    struct bond_angle_spring *bond_angles = m->bond_angles;

    //Bond angle constraints
    #ifdef HAVE_OPENMP
    #pragma omp parallel for shared(bond_angles)
    #endif
    for(size_t i=0; i < m->num_bond_angles; i++){
        struct bond_angle_spring *s = &bond_angles[i];

        if(s->a1->fixed && s->a2->fixed && s->a3->fixed)
            continue;

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
            if(m->debug)
                debug_angle(m, s);
        }
    }
}

void apply_drag_force(struct model *m){
    struct vector tmp;
    //If we're not using the fancy drag force, apply the drag force now.
    if(!m->shield_drag){
        for(size_t i=0; i < m->num_atoms; i++){
            vector_copy_to(&tmp, &m->atoms[i].velocity);
            vmul_by(&tmp, m->drag_coefficient);
            vadd_to(&m->atoms[i].force, &tmp);
        }
    }
}

double model_energy(struct model *m){
    double energy = 0;

    for(size_t i = 0; i < m->num_linear_springs; i++)
        if(linear_spring_synthesised(&m->linear_springs[i]))
            if(linear_spring_active(&m->linear_springs[i]))
                energy += linear_spring_energy(&m->linear_springs[i]);

    for(size_t i = 0; i < m->num_bond_angles; i++)
        if(bond_angle_synthesised(&m->bond_angles[i]))
            energy += bond_angle_energy(&m->bond_angles[i]);

    for(size_t i = 0; i < m->num_torsion_springs; i++)
        if(torsion_spring_synthesised(&m->torsion_springs[i]))
            energy += torsion_spring_energy(&m->torsion_springs[i]);

    return energy;
}

double constraint_energy(struct constraint *c, struct model *m){
    //Model this as a quadratic potential with a high constant
    struct atom *a = &m->atoms[c->a];
    struct atom *b = &m->atoms[c->b];

    struct vector displacement;
    vsub(&displacement, &a->position, &b->position);
    double distance = vmag(&displacement);

    double k = 1;
    double dr = distance - c->distance;
    return k * dr*dr;
}

void constraint_force(struct constraint *c, struct model *m){
    //Model this as a quadratic potential with a high constant
    struct atom *a = &m->atoms[c->a];
    struct atom *b = &m->atoms[c->b];

    struct vector displacement;
    vsub(&displacement, &a->position, &b->position);
    double distance = vmag(&displacement);

    double k = 1;
    double dr = distance - c->distance;

    struct vector force_a, force_b;
    vmul(&force_a, &displacement, -dr * k);
    vmul(&force_b, &force_a, -1);

    vadd_to(&a->force, &force_a);
    vadd_to(&b->force, &force_b);
}

bool constraint_is_synthesised(struct constraint *c, struct model *m){
    struct atom *a = &m->atoms[c->a];
    struct atom *b = &m->atoms[c->b];
    return a->synthesised && b->synthesised;
}

static int m_i = 0;
void model_minim(struct model *m){
    double precision = 0.01;
    double movement = DBL_MAX;

    //We are going to do an inexact line search.

    //Constant in Armijo–Goldstein condition
    const double c = 0.5;

    //Multiplicative factor when reducing step size
    const double tau = 0.5;

    //Energy of system after proposed step
    double curr_energy;
    //Total distance moved by all atoms
    double moved;

    int n = 0;
    do {
        //Just bail out if we try more than a hundred moves
        if(++n > 100)
            break;

        //Get the direction to move in. We are just going to use the negative
        //gradient for a simple steepest descent algorithm.
        struct vector p[m->num_atoms];
        //We will also store the initial position.
        struct vector r[m->num_atoms];

        //Get the initial energy so we can compare it to the energy after any
        //proposed movements to see if the Wolfe conditions are satisfied. More
        //specifically, because we are using the backtracking line search, we check
        //to see if the Armijo–Goldstein condition is satisfied.
        double init_energy = model_energy(m);

        //Get the force. This is equivalent to the negative gradient. Store it in p.
        //We also need to get \f[ \vec{p}^{T} \grad f(\vec{x}) \f]. We are just
        //using the negative gradient, so this is just the dot product \f[ \vec{p}
        //\cdot \vec{p}. \f]

        double total_move = 0;
        double pdot = 0;
        model_accumulate_forces(m);

        //Manually add force due to constraints
        for(size_t i=0; i < m->num_constraints; i++){
            if(constraint_is_synthesised(&m->constraints[i], m)){
                constraint_force(&m->constraints[i], m);
                init_energy += constraint_energy(&m->constraints[i], m);
            }
        }

        for(size_t i=0; i < m->num_atoms; i++){
            vector_copy_to(&p[i], &m->atoms[i].force);
            vector_copy_to(&r[i], &m->atoms[i].position);
            double p_sq = vmag_sq(&p[i]);
            pdot += p_sq;
            total_move += sqrt(p_sq);
        }
        //Get an initial step size
        double step_size = 1;

        while(true){
            model_move_along_vector(m, step_size, r, p);
            curr_energy = model_energy(m);
            if(init_energy - curr_energy < - step_size * c * pdot)
                step_size *= tau;
            else
                break;
        }
        moved = total_move * step_size;
    }while(moved > precision);
}

//Move the model along the direction vector p (with step size alpha), starting at position r.
void model_move_along_vector(struct model *m, double alpha,
        struct vector *r, struct vector *p){

    for(size_t i=0; i < m->num_atoms; i++){
        struct vector dr;
        vector_copy_to(&m->atoms[i].position, &r[i]);
        vector_copy_to(&dr, &p[i]);
        vmul_by(&dr, alpha);
        vadd_to(&m->atoms[i].position, &dr);
    }
}

//Convenience function for profiling to avoid typing the ifdef out
void profile(struct model *m, const char *msg){
    #ifdef HAVE_CLOCK_GETTIME
    if(m->profiler)
        profile_end(m->profiler, "%g\t%s\t%lld\n", m->time, msg);
    #endif
}

void model_build_bond_map(struct model *m){
    //Build a proper 2d array (because we might have a reduced num_atoms to
    //throw off any row/col calculations)
    m->bond_map = malloc(sizeof(bool*) * m->num_atoms);
    for(int i=0; i < m->num_atoms; i++){
        m->bond_map[i] = malloc(sizeof(bool) * m->num_atoms);
        for(int j=0; j < m->num_atoms; j++){
            m->bond_map[i][j] = false;
        }
    }

    //Now add bonds for constraints
    for(int i=0; i < m->num_constraints; i++){
        struct constraint *c = m->constraints + i;
        m->bond_map[c->a][c->b] = true;
        m->bond_map[c->b][c->a] = true;
    }
}

bool model_is_bonded(struct model *m, int i, int j){
    if(!m->bond_map)
        return false;
    return m->bond_map[i][j];
}
