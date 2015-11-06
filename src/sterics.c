#include <stdlib.h>
#include <math.h>
#include "sterics.h"

#define cube(x)   ((x) * (x) * (x))
#define square(x) ((x) * (x))

unsigned long long called = 0;
static int min(int a, int b);
static int max(int a, int b);

//only needed because the generic callback passed to steric_grid_foreach_nearby
//ought to take an extra parameter. This just calls steric_force.
static bool steric_force_lambda(struct atom *a, struct atom *b, void *na);

static bool is_kick_good(struct atom *a, struct atom *b, void *data);
struct is_kick_good_args {
    bool *is_good;
    struct vector *kick_point;
};
static bool should_apply_drag(struct atom *a, struct atom *b, void *data);
static size_t coords_to_index(struct steric_grid *grid, int coords[3]);

struct steric_grid *steric_grid_alloc(double cell_size){
    struct steric_grid *sg = malloc(sizeof(struct steric_grid));
    steric_grid_init(sg, cell_size);
    return sg;
}

void steric_grid_init(struct steric_grid *grid, double cell_size){
    vector_zero(&grid->size);
    grid->cell_size = cell_size;
    grid->atom_grid = NULL;
}

void steric_grid_free(struct steric_grid *grid){
}

/** Generate a lattice with all atoms arranged on it.
 *
 * This will be used for calculating short-range forces like steric effects.
 * The grid will be split into cells each containing a list of atoms falling
 * into that cell. The size of the cells is fixed and chosen to correspond to
 * the steric cut-off.
 *
 * To update the grid, we first scan through the atoms and work out the minimum
 * and maximum values in all dimensions. This allows us to determine how many
 * cells we need in each dimension. If the memory allocated to this steric grid
 * is insufficient, we realloc it to expand.
 */
void steric_grid_update(struct steric_grid *grid, struct model *model){
    vector_fill(&grid->min, INFINITY, INFINITY, INFINITY);
    vector_fill(&grid->max, -INFINITY, -INFINITY, -INFINITY);

    //Find maximum extent of grid in all directions
    for(size_t i=0; i < model->num_atoms; i++){
        vmin_elems(&grid->min, &model->atoms[i].position);
        vmax_elems(&grid->max, &model->atoms[i].position);
    }
    //Bail if no atoms have been synthed
    if(model->num_atoms == 0)
        return;

    //It's best to add a small buffer to the edges here to avoid weird
    //singularities
    for(size_t i=0; i < N; i++){
        grid->max.c[i] += 0.01;
        grid->min.c[i] -= 0.01;
    }

    //Find the number of cells required in each dimension and total number
    size_t num_old_cells = 1;
    size_t num_new_cells = 1;
    for(size_t i=0; i < N; i++){
        num_old_cells *= grid->size.c[i];

        double extent = (grid->max.c[i] - grid->min.c[i]);
        grid->size.c[i] = ceil(extent / grid->cell_size);
        num_new_cells *= grid->size.c[i];
    }
    //Free old buckets
    for(size_t i=0; i < num_old_cells; i++)
        steric_grid_free_atoms(grid, i);

    //Realloc if necessary
    grid->ncells = num_new_cells;
    if(num_new_cells > num_old_cells)
        grid->atom_grid = realloc(grid->atom_grid,
                sizeof(*grid->atom_grid) * num_new_cells);

    //Zero out each grid cell
    for(size_t i=0; i < num_new_cells; i++)
        grid->atom_grid[i] = NULL;

    for(size_t i=0; i < model->num_atoms; i++){
        struct atom *a = &model->atoms[i];
        steric_grid_add_atom(grid, a);
    }
}

//Free linked list of atoms at cell
void steric_grid_free_atoms(struct steric_grid *grid, size_t index){
    struct atom_in_cell *tmp;
    while(grid->atom_grid[index]){
        tmp = grid->atom_grid[index];
        grid->atom_grid[index] = grid->atom_grid[index]->next;
        free(tmp);
    }
}

//Add atom to grid. We alloc a new atom_in_cell struct and add it to the head
//of the linked list.
void steric_grid_add_atom(struct steric_grid *grid, struct atom *a){
    //Find coords

    size_t grid_index = steric_grid_index(grid, a);
    struct atom_in_cell *container = malloc(sizeof(struct atom_in_cell));
    container->a = a;
    container->next = grid->atom_grid[grid_index];
    grid->atom_grid[grid_index] = container;
}

size_t steric_grid_index(struct steric_grid *grid, struct atom *a){
    int coords[N];
    steric_grid_coords(grid, a, coords);
    return coords_to_index(grid, coords);
}

size_t coords_to_index(struct steric_grid *grid, int coords[3]){
    size_t index = coords[0]
        + grid->size.c[0] * coords[1]
        + grid->size.c[0] * grid->size.c[1] * coords[2];
    return index;
}


void steric_grid_coords(struct steric_grid *grid, struct atom *a, int *dst){
    /*
     * Consider an atom a on a line:
     *
     *             a
     * |----|----|----|----|----|----|
     * min       <--d->             max
     *
     * The grid cell of a is 2, given by floor( (x(a) - min) / d).
     */
    for(size_t i=0; i < N; i++)
        dst[i] = (int)((a->position.c[i] - grid->min.c[i]) / grid->cell_size);
}

void steric_grid_forces(struct steric_grid *grid, struct model *model){
    for(size_t i=0; i < model->num_atoms; i++){
        //Get the list of nearby atoms; that is, all atoms in this cell or
        //the surrounding cells
        struct atom *a = &model->atoms[i]; //grid->atom_grid[i][j];
        if(a->fixed)
            continue;
        steric_grid_foreach_nearby(grid, a, steric_force_lambda, NULL);
    }
}

void steric_grid_foreach_nearby(
        struct steric_grid *grid, struct atom *a,
        bool (*lambda)(struct atom *a, struct atom *b, void *data),
        void *data){
    int r[3];
    steric_grid_coords(grid, a, r);

    for(int i=max(0, r[0]-1); i <= min(grid->size.c[0]-1, r[0]+1); i++){
        for(int j=max(0, r[1]-1); j <= min(grid->size.c[1]-1, r[1]+1); j++){
            for(int k=max(0, r[2]-1); k <= min(grid->size.c[2]-1, r[2]+1); k++){
                int coords[3] = {i, j, k};
                int index = coords_to_index(grid, coords);

                struct atom_in_cell *c;
                for(c = grid->atom_grid[index]; c; c = c->next){
                    if(c->a == a)
                        continue;

                    bool cont = (*lambda)(a, c->a, data);
                    if(!cont)
                        return;
                }
            }
        }
    }
}

static bool steric_force_lambda(struct atom *a, struct atom *b, void *na){
    steric_force(a, b);
    return true;
}

void steric_force(struct atom *a, struct atom *b){
    called++;
    struct vector displacement;
    vsub(&displacement, &a->position, &b->position);
    double dist = vmag(&displacement);
    if(dist < a->radius + b->radius){
        double excess = (dist - (a->radius + b->radius));
        vmul_by(&displacement,
                STERIC_FORCE_CONSTANT
                * excess / dist);


        //Only apply in one direction; because of the way we are iterating over
        //the atoms, this function will be called with the other atom as a
        vmul_by(&displacement, -1);
        vadd_to(&a->force, &displacement);
    }
}

#define POLAR_KICK_PROB 0.0001
#define KICK_PROB 0.0003
#define WATER_RADIUS 1.4
#define DRAG_SHIELDING_DISTANCE 10.14
#define COS_DRAG_BLOCK_ANGLE 0.80901699437494745
#define KICK_VELOCITY 0.08

void water_force(struct model *m, struct steric_grid *grid){
    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a = &m->atoms[i];
        if(a->fixed)
            continue;

        double sf_area = 4*M_PI*a->radius*a->radius;

        double kick_prob = sf_area * (POLAR_KICK_PROB + (a->hydrophobicity
                    * (KICK_PROB - POLAR_KICK_PROB)));

        bool do_kick = kick_prob * RAND_MAX * m->timestep < rand();

        struct vector kick, kick_point;
        vector_rand(&kick, 0, M_PI);

        vector_copy_to(&kick_point, &kick);
        vmul_by(&kick_point, a->radius);
        vadd_to(&kick_point, &a->position);

        if(do_kick){
            bool good = true;
            struct is_kick_good_args args = {&good, &kick_point};
            steric_grid_foreach_nearby(grid, a, is_kick_good, &args);
            if(good){
                vmul_by(&kick, KICK_VELOCITY);
                vadd_to(&a->force, &kick);
            }
        }
    }
}

void drag_force(struct model *m, struct steric_grid *grid){
    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a = &m->atoms[i];
        if(a->fixed)
            continue;

        struct vector drag;
        vector_copy_to(&drag, &a->velocity);

        if(vmag(&drag) > 0){
            bool apply_drag = true;
            steric_grid_foreach_nearby(grid, a, should_apply_drag, &apply_drag);
            if(apply_drag){
                vmul_by(&drag, m->drag_coefficient);
                vadd_to(&a->force, &drag);
            }
        }
    }
}

bool is_kick_good(struct atom *a, struct atom *b, void *data){
    struct is_kick_good_args *args = (struct is_kick_good_args *) data;

    struct vector displ;
    vsub(&displ, args->kick_point, &b->position);
    if(vmag(&displ) < b->radius + WATER_RADIUS){
        (*args->is_good) = false;
        return false;
    }
    return true;
}

bool should_apply_drag(struct atom *a, struct atom *b, void *data){
    bool *apply = (bool*) data;

    struct vector displ;
    vsub(&displ, &b->position, &a->position);
    double dist = vmag(&displ);

    if(dist < DRAG_SHIELDING_DISTANCE){
        double dot = vdot(&displ, &a->velocity);
        if(dot / (dist * vmag(&a->velocity)) > COS_DRAG_BLOCK_ANGLE){
            (*apply) = false;
            return false;
        }
    }
    return true;
}

int min(int a, int b){
    return (a < b) ? a : b;
}

int max(int a, int b){
    return (a > b) ? a : b;
}
