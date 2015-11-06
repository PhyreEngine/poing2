#include <stdlib.h>
#include "leapfrog.h"
#include "rk4.h"
#include "vector.h"

//The leapfrog integrator requires the velocity to be a half-step out of phase
//with the position (hence "leapfrog"). To initialise it, we use the RK4
//integrator with a half timestep.
void leapfrog_init(struct model *model){
    double dt = model->timestep;

    size_t num_atoms = 0;
    for(size_t i=0; i < model->num_residues; i++)
        num_atoms += model->residues[i].num_atoms;

    //Store original positions; we only want to increase the velocity
    struct vector positions[num_atoms];

    size_t k = 0;
    for(size_t i=0; i < model->num_residues; i++)
        for(size_t j=0; j < model->residues[i].num_atoms; j++)
            vector_copy_to(&positions[k++], &model->residues[i].atoms[j].position);

    model_accumulate_forces(model);
    model->timestep = dt / 2;
    rk4_push(model);
    model->timestep = dt;

    k = 0;
    for(size_t i=0; i < model->num_residues; i++)
        for(size_t j=0; j < model->residues[i].num_atoms; j++)
            vector_copy_to(&model->residues[i].atoms[j].position, &positions[k++]);
}

void leapfrog_push(struct model *model){
    //Temp vectors for storing change in position and velocity.
    //We use two separate vectors for clarity.
    struct vector dr, dv;

    double dt = model->timestep;
    model->time += dt;

    //Build a flat array of atoms, just to make life a bit simpler
    size_t num_atoms = 0;
    for(size_t i=0; i < model->num_residues; i++)
        num_atoms += model->residues[i].num_atoms;

    struct atom *atoms[num_atoms];
    size_t k = 0;
    for(size_t i=0; i < model->num_residues; i++)
        for(size_t j=0; j < model->residues[i].num_atoms; j++)
            atoms[k++] = &model->residues[i].atoms[j];

    //Calculate forces
    model_accumulate_forces(model);

    //Increase position and velocity
    for(size_t i=0; i < num_atoms; i++){
        vmul(&dr, &atoms[i]->velocity, dt);
        vadd_to(&atoms[i]->position, &dr);

        vmul(&dv, &atoms[i]->force, dt / atoms[i]->mass);
        vadd_to(&atoms[i]->velocity, &dv);
    }
}
