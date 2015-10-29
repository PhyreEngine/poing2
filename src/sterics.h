#ifndef STERICS_H_
#define STERICS_H_

#define GRID_BUFFER 0.01
#define MAX_STERIC_DISTANCE 5.0
#define STERIC_FORCE_CONSTANT 1.025

#include "residue.h"
#include "model.h"
#include "vector.h"

//The grid is essentially a hash table, and we are nearly guaranteed to have
//collisions (it's a feature). We use a linked list to store the atoms in each
//cell.
struct atom_in_cell {
    struct atom *a;
    struct atom_in_cell *next;
};

struct steric_grid {
    struct vector min, max;
    ///Length, width and height
    struct vector size;
    double cell_size;
    struct atom_in_cell **atom_grid;

    size_t ncells;
};

struct steric_grid *steric_grid_alloc(double cell_size);
void steric_grid_init(struct steric_grid *grid, double cell_size);
void steric_grid_free(struct steric_grid *grid);
void steric_grid_update(struct steric_grid *grid, struct model *model);
void steric_grid_forces(struct steric_grid *grid, struct model *model);
size_t steric_grid_index(struct steric_grid *grid, struct atom *a);
void steric_grid_coords(struct steric_grid *grid, struct atom *a, int *dst);
void steric_force(struct atom *a, struct atom *b);
void steric_grid_free_atoms(struct steric_grid *grid, size_t index);

void steric_grid_add_atom(struct steric_grid *grid, struct atom *a);

void steric_grid_foreach_nearby(
        struct steric_grid *grid, struct atom *a,
        bool (*lambda)(struct atom *a, struct atom *b, void *data),
        void *data);

void water_force(struct model *m, struct steric_grid *grid);
void drag_force(struct model *m, struct steric_grid *grid);

#endif /* STERICS_H_ */

