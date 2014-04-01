#ifndef MODEL_H_
#define MODEL_H_

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"

/**
 * Represents a model of a protein, with residues and springs.
 */

struct model {
    ///Number of linear springs.
    size_t num_linear_springs;
    ///Number of torsion springs.
    size_t num_torsion_springs;
    ///Number of residues
    size_t num_residues;

    ///Residues
    struct residue *residues;
    ///Linear springs
    struct linear_spring *linear_springs;
    ///Torsion springs
    struct torsion_spring *torsion_springs;

    ///Current time
    double time;
    ///Timestep
    double timestep;
    ///Time between residues being synthesised
    double synth_time;
    ///Drag coefficient
    double drag_coefficient;
};

struct model *model_alloc();
void model_free(struct model *m);

void model_accumulate_forces(struct model *m);
int model_pdb(FILE *out, const struct model *m, bool conect);
void model_synth(struct model *state, const struct model *m);


#endif /* MODEL_H_ */

