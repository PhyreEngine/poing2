#include <stdlib.h>
#include "leapfrog.h"
#include "rk4.h"
#include "vector.h"
#include "model.h"
#include "residue.h"

//The leapfrog integrator requires the velocity to be a half-step out of phase
//with the position (hence "leapfrog"). To initialise it, we use the RK4
//integrator with a half timestep.
void leapfrog_init(struct model *model){
    double dt = model->timestep;

    //Store original positions; we only want to increase the velocity
    struct vector positions[model->num_atoms];

    for(size_t i=0; i < model->num_atoms; i++)
        vector_copy_to(&positions[i], &model->atoms[i].position);

    model_accumulate_forces(model);
    model->timestep = dt / 2;
    rk4_push(model);
    model->timestep = dt;

    for(size_t i=0; i < model->num_atoms; i++)
        vector_copy_to(&model->atoms[i].position, &positions[i]);
}

void leapfrog_push(struct model *model){
    //Temp vectors for storing change in position and velocity.
    //We use two separate vectors for clarity.
    struct vector dr, dv;

    double dt = model->timestep;
    model->time += dt;

    //Calculate forces
    model_accumulate_forces(model);

    //Increase position and velocity
    for(size_t i=0; i < model->num_atoms; i++){
        if(model->atoms[i].fixed)
            continue;

        vmul(&dr, &model->atoms[i].velocity, dt);
        vadd_to(&model->atoms[i].position, &dr);

        vmul(&dv, &model->atoms[i].force, dt / model->atoms[i].mass);
        vadd_to(&model->atoms[i].velocity, &dv);
    }
}
