#ifndef MODEL_H_
#define MODEL_H_

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>

struct residue;

#define DEFAULT_MAX_SYNTH_ANGLE 1
struct steric_grid;

#define DEBUG_LINEAR_FIELDS \
    "time " \
    "a.id a.name b.id b.name " \
    "enabled " \
    "equilibrium distance " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z\n"
#define DEBUG_LINEAR_FMT \
    "%f " \
    "%d %s %d %s " \
    "%s " \
    "%f %f " \
    "%f %f %f " \
    "%f %f %f\n"

#define DEBUG_ANGLE_FIELDS \
    "time " \
    "a.id a.name b.id b.name c.id c.name " \
    "enabled " \
    "equilibrium angle " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z " \
    "force.c.x force.c.y force.c.z\n"
#define DEBUG_ANGLE_FMT \
    "%f " \
    "%d %s %d %s %d %s " \
    "%s " \
    "%f %f "\
    "%f %f %f %f %f %f %f %f %f\n"

#define DEBUG_TORSION_FIELDS \
    "time " \
    "a.id a.name b.id b.name c.id c.name d.id d.name " \
    "enabled " \
    "equilibrium angle " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z " \
    "force.c.x force.c.y force.c.z " \
    "force.d.x force.d.y force.d.z\n"
#define DEBUG_TORSION_FMT \
    "%f " \
    "%d %s %d %s %d %s %d %s " \
    "%s " \
    "%f %f "\
    "%f %f %f %f %f %f %f %f %f %f %f %f\n"

struct constraint {
    //Atom indices
    size_t a, b;
    float distance;
};

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
    ///Number of atoms
    size_t num_atoms;
    ///Number of bond angle constraints.
    size_t num_bond_angles;
    ///Number of Ramachandran constraints
    size_t num_rama_constraints;
    ///Number of hard constraints
    size_t num_constraints;

    ///Residues
    struct residue *residues;
    ///Atoms
    struct atom *atoms;
    ///Linear springs
    struct linear_spring *linear_springs;
    ///Torsion springs
    struct torsion_spring *torsion_springs;
    ///Bond angle springs
    struct bond_angle_spring *bond_angles;
    ///Ramachandran constraints
    struct rama_constraint *rama_constraints;
    ///Hard constraints
    struct constraint *constraints;

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

    ///Fix all atoms before L-n if the jitter is low, where L is number of
    //currently synthesised residues and n is this variablep.
    int fix_before;

    ///Record the position at this time step;
    double record_time;

    ///Freeze atoms with jitter below this
    double max_jitter;

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

    ///Optional files to write debugging information to
    struct model_debug *debug;
};

///Different output files for debugging information
struct model_debug {
    ///File to print debug info about linear springs
    FILE *linear;
    ///File to print debug info about bond angle springs
    FILE *angle;
    ///File to print debug info about torsion springs
    FILE *torsion;
    ///File to print debug info about Ramachandran constraints
    FILE *rama;

    ///Print results at this timestep interval
    double interval;
    ///Number of times we have printed
    size_t nprinted;
};

struct model *model_alloc();
void model_free(struct model *m);

void model_accumulate_forces(struct model *m);
int model_pdb(FILE *out, const struct model *m, bool conect);
void model_synth(struct model *state, const struct model *m);

void model_synth_atom(const struct model *m, size_t idx, double max_angle);

#endif /* MODEL_H_ */

