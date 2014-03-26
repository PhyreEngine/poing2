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
    //Begin by zeroing out any existing forces
    for(size_t i=0; i < m->num_residues; i++)
        vector_fill(m->residues[i].force, 0, 0, 0);

    //Then go through all springs and accumulate forces on the residues
    for(size_t i=0; i < m->num_linear_springs; i++){
        struct linear_spring s = m->linear_springs[i];
        vadd_to(s.a->force, linear_spring_force(&s, A));
        vadd_to(s.b->force, linear_spring_force(&s, B));
    }

    for(size_t i=0; i < m->num_torsion_springs; i++){
        struct torsion_spring s = m->torsion_springs[i];
        vadd_to(s.r1->force, torsion_spring_force(&s, R1));
        vadd_to(s.r4->force, torsion_spring_force(&s, R4));
    }
}
