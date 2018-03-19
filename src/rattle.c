#include "vector.h"
#include "model.h"
#include "residue.h"
#include "rattle.h"
#include <math.h>
#include <signal.h>

static size_t maxit = 100;
static double tolerance = 1e-4;

static int ni = 0;

void rattle_push(struct model *m){
    rattle_unconstrained_push(m);
    model_accumulate_forces(m);
    rattle_move(m);
    m->time += m->timestep;
}

void rattle_unconstrained_push(struct model *m){
    ni++;
    bool moving[m->num_atoms];
    bool moved[m->num_atoms];

    //We will need to store the unconstrained position of each atom after the
    //initial push.
    struct vector uncons[m->num_atoms];

    //Do the initial verlet push, storing the positions in "ucons"
    for(size_t a=0; a < m->num_atoms; a++){
        if(m->atoms[a].fixed){
            vector_copy_to(&uncons[a], &m->atoms[a].position);
            vector_zero(&m->atoms[a].velocity);
            continue;
        }

        //Get acceleration
        struct vector accel;
        vdiv(&accel, &m->atoms[a].force, m->atoms[a].mass);

        moving[a] = false;
        moved[a] = true;

        //Get unconstrained position by a velocity Verlet push
        for(size_t i=0; i < N; i++){
            uncons[a].c[i] = m->atoms[a].position.c[i]
                + m->timestep * m->atoms[a].velocity.c[i]
                + m->timestep * m->timestep / 2 * accel.c[i];
            //Also push the velocity a half step
            m->atoms[a].velocity.c[i] =
                m->atoms[a].velocity.c[i] + m->timestep / 2 * accel.c[i];
        }
    }


    //Begin iterating to solve the constraints
    bool done = false;
    for(size_t nit = 0; !done && nit < maxit; nit++){
        for(size_t i=0; i < m->num_constraints; i++){
            //Set to false if anything is moved.
            done = true;

            struct atom *a = &m->atoms[m->constraints[i].a];
            struct atom *b = &m->atoms[m->constraints[i].b];
            if(!a->synthesised || !b->synthesised)
                continue;
            if(!moved[m->constraints[i].a] && !moved[m->constraints[i].b])
                continue;

            //Get displacement vector between unconsrained positions
            struct vector p;
            vsub(&p, &uncons[m->constraints[i].a], &uncons[m->constraints[i].b]);

            //Do we need to apply this constaint?
            float dist = m->constraints[i].distance;
            float diffsq = dist * dist - vmag_sq(&p);
            if(fabs(diffsq) > tolerance * 2){
                //This is not our last iteration
                done = false;

                //Get displacement vector between unmoved atoms
                struct vector r;
                vsub(&r, &a->position, &b->position);

                //XXX: The original Allen and Tildsey code has a bail out here
                //if a certain tolerance is not met.

                //Get correction factor g_ab
                float reduced_mass = 1.0 / a->mass + 1.0 / b->mass;
                float gab = diffsq / (2.0 * reduced_mass * vdot(&r, &p));

                //Get correction term
                struct vector delta;
                vmul(&delta, &r, gab);

                //Update unconstrained positions
                for(size_t j=0; j<N; j++){
                    if(!m->atoms[m->constraints[i].a].fixed){
                        uncons[m->constraints[i].a].c[j] += 1.0 / a->mass * delta.c[j];
                        a->velocity.c[j] += 1.0 / a->mass * delta.c[j] / m->timestep;
                    }
                    if(!m->atoms[m->constraints[i].b].fixed){
                        uncons[m->constraints[i].b].c[j] -= 1.0 / b->mass * delta.c[j];
                        b->velocity.c[j] -= 1.0 / b->mass * delta.c[j] / m->timestep;
                    }
                }
            }
        }
        for(size_t i=0; i < m->num_atoms; i++){
            moved[i] = moving[i];
            moving[i] = false;
        }
    }
    if(!done)
        fprintf(stderr, "Warning: Maximum iterations exceeded at line %d of file %s\n", __LINE__, __FILE__);

    //Copy the new positions to the atoms
    for(size_t i=0; i < m->num_atoms; i++){
        if(!m->atoms[i].fixed)
            vector_copy_to(&m->atoms[i].position, &uncons[i]);
    }
}

size_t ncalled = 0;

//Call this after calculating new forces
void rattle_move(struct model *m){ ncalled++;
    bool moving[m->num_atoms];
    bool moved[m->num_atoms];

    if(getenv("DEBUG"))
        if(ncalled == 235001)
            raise(SIGINT);

    //Do the second verlet push
    for(size_t a=0; a < m->num_atoms; a++){
        if(m->atoms[a].fixed)
            continue;

        //Update velocity using acceleration
        for(size_t i=0; i < N; i++){
            m->atoms[a].velocity.c[i] =
                m->atoms[a].velocity.c[i]
                + (m->timestep / 2) * m->atoms[a].force.c[i] / m->atoms[a].mass;
        }

        moving[a] = false;
        moved[a] = true;
    }

    //Begin iterating to converge on velocity
    bool done = false;
    for(size_t nit = 0; nit < maxit && !done; nit++){
        done = true;
        for(size_t i=0; i < m->num_constraints; i++){
            struct atom *a = &m->atoms[m->constraints[i].a];
            struct atom *b = &m->atoms[m->constraints[i].b];
            if(!a->synthesised || !b->synthesised)
                continue;
            if(!moved[m->constraints[i].a] && !moved[m->constraints[i].b])
                continue;

            //Constraint distance squared
            float dsq = m->constraints[i].distance * m->constraints[i].distance;

            //Get velocity and position delta
            struct vector v_ab, r_ab;
            vsub(&v_ab, &a->velocity, &b->velocity);
            vsub(&r_ab, &a->position, &b->position);

            //Dot product
            float rv = vdot(&v_ab, &r_ab);

            //Reciprocal masses
            float rma = 1.0 / a->mass;
            float rmb = 1.0 / b->mass;

            //Correction term
            float gab = -rv / ((rma + rmb) * dsq);

            if(fabs(gab) > tolerance){
                //Get position correction
                struct vector r_delta;
                vmul(&r_delta, &r_ab, gab);

                //Update velocity vectors
                struct vector delta_v;
                vmul(&delta_v, &r_delta, rma);
                if(!a->fixed)
                    vadd_to(&a->velocity, &delta_v);
                vmul(&delta_v, &r_delta, -rmb);
                if(!b->fixed)
                    vadd_to(&b->velocity, &delta_v);

                done = false;
                moving[m->constraints[i].a] = true;
                moving[m->constraints[i].b] = true;
            }
        }

        for(size_t i=0; i < m->num_atoms; i++){
            moved[i] = moving[i];
            moving[i] = false;
        }
    }
    if(!done)
        fprintf(stderr, "Warning: Maximum iterations exceeded at line %d of file %s\n", __LINE__, __FILE__);
}
