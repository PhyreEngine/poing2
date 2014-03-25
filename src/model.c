#include <stdlib.h>
#include "model.h"

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
