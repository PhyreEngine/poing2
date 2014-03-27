#include <stdlib.h>
#include <stdio.h>
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

    //Begin by zeroing out any existing forces
    for(size_t i=0; i < m->num_residues; i++)
        vector_zero(&m->residues[i].force);

    //Then go through all springs and accumulate forces on the residues
    for(size_t i=0; i < m->num_linear_springs; i++){
        struct linear_spring s = m->linear_springs[i];

        linear_spring_force(&force, &s, A);
        vadd_to(&s.a->force, &force);

        linear_spring_force(&force, &s, B);
        vadd_to(&s.b->force, &force);
    }

    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring s = m->torsion_springs[i];

        torsion_spring_force(&force, &s, R1);
        vadd_to(&s.r1->force, &force);

        torsion_spring_force(&force, &s, R4);
        vadd_to(&s.r4->force, &force);
    }

    for(size_t i=0; i < m->num_residues; i++){
        vector_copy_to(&tmp, &m->residues[i].velocity);
        vmul_by(&tmp, -1);
        vadd_to(&m->residues[i].force, &tmp);
    }
}


const char *fmt = "ATOM  %5d %4s %-3s  %4d%1s   %8.3f%8.3f%8.3f\n";
char * model_pdb(const struct model *m){
    char *buffer = malloc(m->num_residues * 80);
    buffer[0] = '\0';

    char line[80];
    for(size_t i=0; i < m->num_residues; i++){
        sprintf(line, fmt, i+1, " CA ", 
                m->residues[i].aa->threeletter,
                i+1, " ",
                m->residues[i].position.c[0],
                m->residues[i].position.c[1],
                m->residues[i].position.c[2]);
        strcat(buffer, line);
    }
    return buffer;
}
