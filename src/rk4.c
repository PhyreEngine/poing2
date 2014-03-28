#include <stdlib.h>
#include "rk4.h"
#include "vector.h"

void rk4_push(struct model *model){
    double dt = model->timestep;
    model->time += dt;

    struct vector orig_pos[model->num_residues];
    struct vector k1[model->num_residues];
    struct vector k2[model->num_residues];
    struct vector k3[model->num_residues];
    struct vector k4[model->num_residues];

    struct residue * residues = model->residues;

    model_accumulate_forces(model);
    #pragma omp parallel for shared(residues, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_residues; i++){
        vector_copy_to(&orig_pos[i], &residues[i].position);
        vector_copy_to(&k1[i], &residues[i].force);

        vector_copy_to(&residues[i].position, &k1[i]);
        vmul_by(&residues[i].position, dt/2);
        vadd_to(&residues[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(residues, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_residues; i++){
        vector_copy_to(&k2[i], &residues[i].force);

        vector_copy_to(&residues[i].position, &k2[i]);
        vmul_by(&residues[i].position, dt/2);
        vadd_to(&residues[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(residues, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_residues; i++){
        vector_copy_to(&k3[i], &residues[i].force);

        vector_copy_to(&residues[i].position, &k3[i]);
        vmul_by(&residues[i].position, dt);
        vadd_to(&residues[i].position, &orig_pos[i]);
    }

    model_accumulate_forces(model);
    #pragma omp parallel for shared(residues, orig_pos, k1, k2, k3, k4)
    for(size_t i=0; i < model->num_residues; i++){
        vector_copy_to(&k4[i], &residues[i].force);

        vmul_by(&k1[i], 1*dt/6);
        vmul_by(&k2[i], 2*dt/6);
        vmul_by(&k3[i], 2*dt/6);
        vmul_by(&k4[i], 1*dt/6);
        vadd_to(&residues[i].velocity, &k1[i]);
        vadd_to(&residues[i].velocity, &k2[i]);
        vadd_to(&residues[i].velocity, &k3[i]);
        vadd_to(&residues[i].velocity, &k4[i]);

        vmul(&residues[i].position, &residues[i].velocity, dt);
        vadd_to(&residues[i].position, &orig_pos[i]);
    }
}
