#include <stdlib.h>
#include "rk4.h"
#include "vector.h"

void rk4_push(struct model *model, double dt){
    vector orig_pos[model->num_residues];
    vector k1[model->num_residues];
    vector k2[model->num_residues];
    vector k3[model->num_residues];
    vector k4[model->num_residues];

    model_accumulate_forces(model);
    for(size_t i=0; i < model->num_residues; i++){
        orig_pos[i] = vector_copy(model->residues[i].position);
        k1[i] = vector_copy(model->residues[i].force);
        vector_fill(model->residues[i].position,
            orig_pos[i][0] + k1[i][0] * dt/2,
            orig_pos[i][1] + k1[i][1] * dt/2,
            orig_pos[i][2] + k1[i][2] * dt/2);
    }

    model_accumulate_forces(model);
    for(size_t i=0; i < model->num_residues; i++){
        k2[i] = vector_copy(model->residues[i].force);
        vector_fill(model->residues[i].position,
            orig_pos[i][0] + k2[i][0] * dt/2,
            orig_pos[i][1] + k2[i][1] * dt/2,
            orig_pos[i][2] + k2[i][2] * dt/2);
    }

    model_accumulate_forces(model);
    for(size_t i=0; i < model->num_residues; i++){
        k3[i] = vector_copy(model->residues[i].force);
        vector_fill(model->residues[i].position,
            orig_pos[i][0] + k3[i][0] * dt,
            orig_pos[i][1] + k3[i][1] * dt,
            orig_pos[i][2] + k3[i][2] * dt);
    }

    model_accumulate_forces(model);
    for(size_t i=0; i < model->num_residues; i++){
        k4[i] = vector_copy(model->residues[i].force);
        vmul_by(k1[i], 1*dt/6);
        vmul_by(k2[i], 2*dt/6);
        vmul_by(k3[i], 2*dt/6);
        vmul_by(k4[i], 1*dt/6);
        vadd_to(model->residues[i].velocity, k1[i]);
        vadd_to(model->residues[i].velocity, k2[i]);
        vadd_to(model->residues[i].velocity, k3[i]);
        vadd_to(model->residues[i].velocity, k4[i]);

        model->residues[i].position[0] =
            orig_pos[i][0] + model->residues[i].velocity[0] * dt;
        model->residues[i].position[1] =
            orig_pos[i][1] + model->residues[i].velocity[1] * dt;
        model->residues[i].position[2] =
            orig_pos[i][2] + model->residues[i].velocity[2] * dt;
        free(k1[i]);
        free(k2[i]);
        free(k3[i]);
        free(k4[i]);
        free(orig_pos[i]);
    }
}
