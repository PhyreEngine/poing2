#ifndef STERICS_H_
#define STERICS_H_

#define GRID_BUFFER 0.01
#define MAX_STERIC_DISTANCE 5.0
#define STERIC_FORCE_CONSTANT 0.1

#include "residue.h"
#include "model.h"
#include "vector.h"

struct steric_grid {
    struct vector min, max;
    size_t divisions;
    struct atom ***atom_grid;
    size_t *num_atoms;
};

void steric_grid_init(struct steric_grid *grid, size_t divisions);
void steric_grid_update(struct steric_grid *grid, struct model *model);
void steric_grid_forces(struct steric_grid *grid, struct model *model);
size_t steric_grid_index(struct steric_grid *grid, struct atom *a);
void steric_grid_coords(struct steric_grid *grid, struct atom *a,
        int *x, int *y, int *z);
void steric_force(struct atom *a, struct atom *b);

#endif /* STERICS_H_ */

