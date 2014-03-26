#include <stdlib.h>
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
