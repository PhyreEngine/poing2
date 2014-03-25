#ifndef MODEL_H_
#define MODEL_H_

#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"

/**
 * Represents a model of a protein, with residues and springs.
 */

struct model {
    ///Number of linear springs.
    int num_linear_springs;
    ///Number of torsion springs.
    int num_torsion_springs;
    ///Number of residues
    int num_residues;

    ///Residues
    struct residue *residues;
    ///Linear springs
    struct linear_spring *linear_springs;
    ///Torsion springs
    struct torsion_spring *torsion_springs;
};

struct model *model_alloc();
void model_free(struct model *m);

#endif /* MODEL_H_ */

