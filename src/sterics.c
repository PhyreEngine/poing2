#include <stdlib.h>
#include <math.h>
#include "sterics.h"

static void apply_forces_nearby(struct steric_grid *grid, struct atom *a);
static int min(int a, int b);
static int max(int a, int b);

void steric_grid_init(struct steric_grid *grid, size_t divisions){
    double d3 = divisions * divisions * divisions;
    grid->divisions = divisions;
    grid->atom_grid = malloc(sizeof(struct atom *) * d3);
    grid->num_atoms = malloc(sizeof(struct atom *) * d3);
    for(size_t i=0; i < d3; i++){
        grid->atom_grid[i] = NULL;
        grid->num_atoms[i] = 0;
    }
    vector_fill(&grid->min, 0, 0, 0);
    vector_fill(&grid->max, 0, 0, 0);
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
                    (grid->num_atoms[index] + 1) * sizeof(a));

            grid->atom_grid[index][grid->num_atoms[index]] = a;
            grid->num_atoms[index]++;
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
    size_t d3 = grid->divisions * grid->divisions * grid->divisions;
    for(size_t i=0; i < d3; i++){
        for(size_t j=0; j < grid->num_atoms[i]; j++){
            //Get the list of nearby atoms; that is, all atoms in this cell or
            //the surrounding cells
            struct atom *a = grid->atom_grid[i][j];
            apply_forces_nearby(grid, a);
        }
    }
}

//Ugh, triangle code.
void apply_forces_nearby(struct steric_grid *grid, struct atom *a){
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
                size_t num_atoms    = grid->num_atoms[index];
                for(size_t l=0; l < num_atoms; l++){
                    if(a != atoms[l])
                        steric_force(a, atoms[l]);
                }
            }
        }
    }
}

void steric_force(struct atom *a, struct atom *b){
    struct vector displacement;
    vsub(&displacement, &a->position, &b->position);
    double dist = vmag(&displacement);
    if(dist < a->radius + b->radius){
        vmul_by(&displacement,
                STERIC_FORCE_CONSTANT 
                * (dist - (a->radius + b->radius)) / dist);

        //Only apply in one direction; because of the way we are iterating over
        //the atoms, this function will be called with the other atom as a
        vmul_by(&displacement, -1);
        vadd_to(&a->force, &displacement);
    }
}

int min(int a, int b){
    return (a < b) ? a : b;
}

int max(int a, int b){
    return (a > b) ? a : b;
}
