#include <stdlib.h>
#include "rk4.h"
#include "vector.h"
#include "model.h"
#include "residue.h"

void rk4_push(struct model *model){
    double dt = model->timestep;
    model->time += dt;

    struct vector orig_pos[model->num_atoms];
    struct vector k1[model->num_atoms];
    struct vector k2[model->num_atoms];
    struct vector k3[model->num_atoms];
    struct vector k4[model->num_atoms];

    struct atom *atoms = model->atoms;

    model_accumulate_forces(model);
    #pragma omp parallel for shared(atoms, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_atoms; i++){
        vector_copy_to(&orig_pos[i], &atoms[i].position);
        vector_copy_to(&k1[i], &atoms[i].force);
        vdiv_by(&k1[i], atoms[i].mass);

        vector_copy_to(&atoms[i].position, &k1[i]);
        vmul_by(&atoms[i].position, dt/2);
        vadd_to(&atoms[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(atoms, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_atoms; i++){
        vector_copy_to(&k2[i], &atoms[i].force);
        vdiv_by(&k2[i], atoms[i].mass);

        vector_copy_to(&atoms[i].position, &k2[i]);
        vmul_by(&atoms[i].position, dt/2);
        vadd_to(&atoms[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(atoms, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_atoms; i++){
        vector_copy_to(&k3[i], &atoms[i].force);
        vdiv_by(&k3[i], atoms[i].mass);

        vector_copy_to(&atoms[i].position, &k3[i]);
        vmul_by(&atoms[i].position, dt);
        vadd_to(&atoms[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(atoms, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_atoms; i++){
        vector_copy_to(&k4[i], &atoms[i].force);
        vdiv_by(&k4[i], atoms[i].mass);

        vmul_by(&k1[i], 1*dt/6);
        vmul_by(&k2[i], 2*dt/6);
        vmul_by(&k3[i], 2*dt/6);
        vmul_by(&k4[i], 1*dt/6);
        vadd_to(&atoms[i].velocity, &k1[i]);
        vadd_to(&atoms[i].velocity, &k2[i]);
        vadd_to(&atoms[i].velocity, &k3[i]);
        vadd_to(&atoms[i].velocity, &k4[i]);

        vmul(&atoms[i].position, &atoms[i].velocity, dt);
        vadd_to(&atoms[i].position, &orig_pos[i]);
    }
}
