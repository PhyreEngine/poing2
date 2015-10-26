#include <stdlib.h>
#include <math.h>
#include "sterics.h"

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

struct steric_grid *steric_grid_alloc(size_t divisions){
    struct steric_grid *sg = malloc(sizeof(struct steric_grid));
    steric_grid_init(sg, divisions);
    return sg;
}

void steric_grid_init(struct steric_grid *grid, size_t divisions){
    double d3 = divisions * divisions * divisions;
    grid->divisions = divisions;
    grid->atom_grid = malloc(sizeof(struct atom *) * d3);
    grid->atoms_per_cell = malloc(sizeof(struct atom *) * d3);
    for(size_t i=0; i < d3; i++){
        grid->atom_grid[i] = NULL;
        grid->atoms_per_cell[i] = 0;
    }
    vector_fill(&grid->min, 0, 0, 0);
    vector_fill(&grid->max, 0, 0, 0);
}

void steric_grid_free(struct steric_grid *grid){
    if(!grid)
        return;
    double d3 = grid->divisions * grid->divisions * grid->divisions;
    for(size_t i=0; i < d3; i++)
        free(grid->atom_grid[i]);

    free(grid->atom_grid);
    free(grid->atoms_per_cell);
    free(grid);
}

/**
 * Generate a lattice with all atoms arranged on it.
 *
 * This will be used for calculating short-range forces like steric effects.
 * The grid will be split into cells each containing a list of atoms falling
 * into that cell. The size of the cells is given by the maximum volume taken
 * by the protein divided by the number of divisions specified.
 *
 * To update the grid, we first scan through the atoms and work out the minimum
 * and maximum values in all dimensions. This allows us to calculate the cell
 * size and we then scan through again and assign all atoms to a cell.
 */
void steric_grid_update(struct steric_grid *grid, struct model *model){
    vector_fill(&grid->min, INFINITY, INFINITY, INFINITY);
    vector_fill(&grid->max, -INFINITY, -INFINITY, -INFINITY);

    double d3 = grid->divisions * grid->divisions * grid->divisions;
    for(size_t i=0; i < d3; i++)
        grid->atoms_per_cell[i] = 0;

    for(size_t i=0; i < model->num_residues; i++){
        for(size_t j=0; j < model->residues[i].num_atoms; j++){
            for(size_t k=0; k < N; k++){
                grid->min.c[k] = fmin(
                        grid->min.c[k],
                        model->residues[i].atoms[j].position.c[k]);
                grid->max.c[k] = fmax(
                        grid->max.c[k],
                        model->residues[i].atoms[j].position.c[k]);
            }
        }
    }
    //It's best to add a small buffer to the edges here to avoid weird
    //singularities
    for(size_t i=0; i < N; i++){
        grid->max.c[i] += 0.01;
        grid->min.c[i] -= 0.01;
    }

    for(size_t i=0; i < model->num_residues; i++){
        for(size_t j=0; j < model->residues[i].num_atoms; j++){
            struct atom *a = &model->residues[i].atoms[j];
            int index = steric_grid_index(grid, a);

            grid->atom_grid[index] = realloc(
                    grid->atom_grid[index],
                    (grid->atoms_per_cell[index] + 1) * sizeof(a));

            grid->atom_grid[index][grid->atoms_per_cell[index]] = a;
            grid->atoms_per_cell[index]++;
        }
    }
}

size_t steric_grid_index(struct steric_grid *grid, struct atom *a){
    int x, y, z;
    steric_grid_coords(grid, a, &x, &y, &z);
    double divs = grid->divisions;
    return x * divs * divs + y * divs + z;
}

void steric_grid_coords(struct steric_grid *grid, struct atom *a,
        int *x, int *y, int *z){
    double divs = grid->divisions;
    double length = grid->max.c[0] - grid->min.c[0];
    double width  = grid->max.c[1] - grid->min.c[1];
    double height = grid->max.c[2] - grid->min.c[2];

    *x = (a->position.c[0] - grid->min.c[0]) / length * divs;
    *y = (a->position.c[1] - grid->min.c[1]) / width  * divs;
    *z = (a->position.c[2] - grid->min.c[2]) / height * divs;
}

void steric_grid_forces(struct steric_grid *grid, struct model *model){
    for(size_t i=0; i < model->num_residues; i++){
        for(size_t j=0; j < model->residues[i].num_atoms; j++){
            //Get the list of nearby atoms; that is, all atoms in this cell or
            //the surrounding cells
            struct atom *a = &model->residues[i].atoms[j]; //grid->atom_grid[i][j];
            if(a->fixed)
                continue;
            steric_grid_foreach_nearby(grid, a, steric_force_lambda, NULL);
        }
    }
}

void steric_grid_foreach_nearby(
        struct steric_grid *grid, struct atom *a,
        bool (*lambda)(struct atom *a, struct atom *b, void *data),
        void *data){
    int x, y, z;
    steric_grid_coords(grid, a, &x, &y, &z);

    //Find the number of cells that correspond to a realistic steric distance
    int dx = ceil(MAX_STERIC_DISTANCE * grid->divisions
            / (grid->max.c[0] - grid->min.c[0]));
    int dy = ceil(MAX_STERIC_DISTANCE * grid->divisions
            / (grid->max.c[1] - grid->min.c[1]));
    int dz = ceil(MAX_STERIC_DISTANCE * grid->divisions
            / (grid->max.c[2] - grid->min.c[2]));

    for(int i=max(0, x-dx); i <= min(grid->divisions-1, x+dx); i++){
        for(int j=max(0, y-dy); j <= min(grid->divisions-1, y+dy); j++){
            for(int k=max(0, z-dz); k <= min(grid->divisions-1, z+dz); k++){
                int index = i * grid->divisions * grid->divisions
                    + j * grid->divisions
                    + k;
                struct atom **atoms = grid->atom_grid[index];
                size_t num_atoms    = grid->atoms_per_cell[index];
                for(size_t l=0; l < num_atoms; l++){
                    if(a != atoms[l]){
                        bool cont = (*lambda)(a, atoms[l], data);
                        if(!cont)
                            return;
                    }
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
    for(size_t i=0; i < m->num_residues; i++){
        for(size_t j=0; j < m->residues[i].num_atoms; j++){
            struct atom *a = &m->residues[i].atoms[j];
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
}

void drag_force(struct model *m, struct steric_grid *grid){
    for(size_t i=0; i < m->num_residues; i++){
        for(size_t j=0; j < m->residues[i].num_atoms; j++){
            struct atom *a = &m->residues[i].atoms[j];
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
