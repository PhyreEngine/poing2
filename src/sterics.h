#ifndef STERICS_H_
#define STERICS_H_

#define GRID_BUFFER 0.01
#define MAX_STERIC_DISTANCE 5.0
#define STERIC_FORCE_CONSTANT 0.025

#include "residue.h"
#include "model.h"
#include "vector.h"

struct steric_grid {
    struct vector min, max;
    size_t divisions;
    struct atom ***atom_grid;
    size_t *num_atoms;
};

struct steric_grid *steric_grid_alloc(size_t divisions);
void steric_grid_init(struct steric_grid *grid, size_t divisions);
void steric_grid_free(struct steric_grid *grid);
void steric_grid_update(struct steric_grid *grid, struct model *model);
void steric_grid_forces(struct steric_grid *grid, struct model *model);
size_t steric_grid_index(struct steric_grid *grid, struct atom *a);
void steric_grid_coords(struct steric_grid *grid, struct atom *a,
        int *x, int *y, int *z);
void steric_force(struct atom *a, struct atom *b);

void steric_grid_foreach_nearby(
        struct steric_grid *grid, struct atom *a,
        bool (*lambda)(struct atom *a, struct atom *b, void *data),
        void *data);

void water_force(struct model *m, struct steric_grid *grid);
void drag_force(struct model *m, struct steric_grid *grid);

#endif /* STERICS_H_ */

