#ifndef MODEL_H_
#define MODEL_H_

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"
#include "bond_angle.h"

#define DEFAULT_MAX_SYNTH_ANGLE 45
struct steric_grid;

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
    ///Number of bond angle constraints.
    size_t num_bond_angles;
    ///Number of Ramachandran constraints
    size_t num_rama_constraints;

    ///Residues
    struct residue *residues;
    ///Linear springs
    struct linear_spring *linear_springs;
    ///Torsion springs
    struct torsion_spring *torsion_springs;
    ///Bond angle springs
    struct bond_angle_spring *bond_angles;
    ///Ramachandran constraints
    struct rama_constraint *rama_constraints;

    ///Current time
    double time;
    ///Run until this time
    double until;
    ///Timestep
    double timestep;
    ///Time between residues being synthesised
    double synth_time;
    ///Drag coefficient
    double drag_coefficient;
    ///Whether to use the shielded drag force
    bool shield_drag;

    ///Grid from which steric forces are calculated
    struct steric_grid *steric_grid;
    ///Enable / disable steric grid
    bool use_sterics;
    ///Enable / disable water effect
    bool use_water;

    ///Maximum synthesis angle for new residues
    double max_synth_angle;

    ///Fix residues after allowing them to reach equilibrium
    bool fix;

    /** Fix all atoms for the first 2/3 synth_time. For the first
     * synth_time / 3 use linear springs and not torsion springs.
     * Then, for the next synth_time / 3 use torsional springs
     * only. Then enable them all and unfix atoms.
     *
     * This should get the correct distances but not necessarily
     * the correct handedness. Then we fix the handedness with the
     * linear springs turned off so that the torsion spring does
     * not have to overcome the linear springs. All atoms except
     * the curent one are fixed so that effects don't carry on to
     * the rest of the atoms.*/
    bool threestate;

    ///If this is false, all atoms are just dumped in with the default position
    bool do_synthesis;
};

struct model *model_alloc();
void model_free(struct model *m);

void model_accumulate_forces(struct model *m);
int model_pdb(FILE *out, const struct model *m, bool conect);
void model_synth(struct model *state, const struct model *m);
struct residue * model_push_residue(struct model *m, int id);

#endif /* MODEL_H_ */

