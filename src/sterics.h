#ifndef STERICS_H_
#define STERICS_H_

#define GRID_BUFFER 0.01
#define MAX_STERIC_DISTANCE 5.0
#define STERIC_FORCE_CONSTANT 1.025

#include "residue.h"
#include "model.h"
#include "vector.h"

struct atom_list {
    size_t atom_idx;
    struct atom_list *next;
};

struct steric_grid {
    //Size (in cells) in each direction.
    size_t size_in_cells;

    //Calculated dynamically on each update step.
    struct vector origin;

    //Size (in Angstroms) of each cell.
    double cell_size;

    //We store a list of atoms in each cell.
    struct atom_list **cells;

    //Store an explicit interaction list for each atom. Do this by storing the
    //index of the cell the atom is in.
    size_t *interaction_list;

    //Pool of atom lists to avoid repeated malloc/free cycles
    struct atom_list *list_buf;
    size_t list_buf_idx;
};

void steric_grid_init(struct steric_grid *grid, size_t size_in_cells, double cell_size, size_t num_atoms);
void steric_grid_find_origin(struct steric_grid *g, struct model *m);
void steric_grid_update(struct steric_grid *g, struct model *m);
struct atom_list *steric_grid_interaction_list(struct steric_grid *g, size_t atom_index);
void steric_grid_build_ilists(struct steric_grid *g, struct model *m);
void steric_grid_forces(struct steric_grid *g, struct model *m);

void water_force(struct model *m, struct steric_grid *grid);
void drag_force(struct model *m, struct steric_grid *g);


#endif /* STERICS_H_ */

