#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "sterics.h"

#define cube(x)   ((x) * (x) * (x))
#define square(x) ((x) * (x))

static inline size_t coords(struct steric_grid *grid, size_t x, size_t y, size_t z);
static struct vector buf = {.c = {GRID_BUFFER, GRID_BUFFER, GRID_BUFFER}};

void steric_grid_init(struct steric_grid *grid, size_t size_in_cells, double cell_size, size_t num_atoms){
    grid->size_in_cells = size_in_cells;
    grid->cell_size = cell_size;

    //We need this many cells:
    size_t ncells = cube(size_in_cells);

    //Initialise linked list buffer
    grid->list_buf = malloc(num_atoms * sizeof(*grid->list_buf));
    grid->list_buf_idx = 0;

    grid->interaction_list = malloc(num_atoms * sizeof(*grid->interaction_list));
    grid->cells = malloc(ncells * sizeof(*grid->cells));
    //Initialise everything to have zero atoms
    for(size_t i=0; i < ncells; i++)
        grid->cells[i] = NULL;
    for(size_t i=0; i < num_atoms; i++)
        grid->interaction_list[i] = 0;
}

//Small utility function to convert 3d coordinates into the block array format.
static inline size_t coords(struct steric_grid *grid, size_t x, size_t y, size_t z){
    return x + y * grid->size_in_cells + z * square(grid->size_in_cells);
}

//Function to convert coordinates in Angstroms to grid cell coords
static inline void ang2cell(struct steric_grid *grid, struct vector *v, size_t *x, size_t *y, size_t *z){
    *x = (size_t)((v->c[0] - grid->origin.c[0]) / grid->cell_size);
    *y = (size_t)((v->c[1] - grid->origin.c[1]) / grid->cell_size);
    *z = (size_t)((v->c[2] - grid->origin.c[2]) / grid->cell_size);
}

//Min and max functions
static inline int min(double a, double b){
    return (a < b) ? a : b;
}
static inline int man(double a, double b){
    return (a > b) ? a : b;
}

//Find the origin of the grid by locating the lowest dimension of each atom
void steric_grid_find_origin(struct steric_grid *g, struct model *m){
    vector_fill(&g->origin, DBL_MAX, DBL_MAX, DBL_MAX);
    for(size_t i=0; i < m->num_atoms; i++)
        vmin_elems(&g->origin, &m->atoms[i].position);
    vsub_to(&g->origin, &buf);
}

//Add atom to cell i (in block coords)
void add_atom(struct steric_grid *grid, size_t cell_index, size_t atom_idx){
    //Get a linked list object from the buffer
    struct atom_list *list = grid->list_buf + (grid->list_buf_idx++);

    list->next = grid->cells[cell_index];
    list->atom_idx = atom_idx;
    grid->cells[cell_index] = list;
}

//Set the interaction list for the atom at atom_index (in the model) to the
//list of atoms in the cell at cell_index
static void update_ilist(struct steric_grid *grid, size_t atom_index, size_t cell_index){
    grid->interaction_list[atom_index] = cell_index;
}

void free_linked_list(struct atom_list *list){
    while(list){
        struct atom_list *next = list->next;
        free(list);
        list = next;
    }
}

void steric_grid_update(struct steric_grid *g, struct model *m){
    size_t x, y, z;
    size_t cell_index;

    //Reset buffer
    g->list_buf_idx = 0;

    //Reset each cell
    for(size_t i=0; i < m->num_atoms; i++){
        size_t cell_idx = g->interaction_list[i];
        g->cells[cell_idx] = NULL;
    }

    //Find origin
    steric_grid_find_origin(g, m);

    //Add atoms to cells
    for(size_t i=0; i < m->num_atoms; i++){
        //Angstroms to cell coordinates
        ang2cell(g, &m->atoms[i].position, &x, &y, &z);
        //To block coordinates
        cell_index = coords(g, x, y, z);
        //Add atom to linked list
        add_atom(g, cell_index, i);
    }
}

void steric_grid_build_ilists(struct steric_grid *g, struct model *m){
    size_t x, y, z;
    size_t cell_index;

    for(size_t i=0; i < m->num_atoms; i++){
        //Angstroms to cell coordinates
        ang2cell(g, &m->atoms[i].position, &x, &y, &z);
        //To block coordinates
        cell_index = coords(g, x, y, z);
        //Add atom to linked list
        update_ilist(g, i, cell_index);
    }
}

struct atom_list *steric_grid_interaction_list(struct steric_grid *g, size_t atom_index){
    return g->cells[g->interaction_list[atom_index]];
}

void steric_grid_forces(struct steric_grid *g, struct model *m){
    struct vector displacement;
    struct atom_list *l = NULL;
    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a = m->atoms + i;

        for(l = g->cells[g->interaction_list[i]]; l; l = l->next){
            struct atom *b = m->atoms + l->atom_idx;
            if(a == b)
                continue;

            vsub(&displacement, &b->position, &a->position);
            double dist = vmag(&displacement);
            if(!model_is_bonded(m, i, l->atom_idx) && dist < a-> radius + b->radius){
                //Find the distance by which the constraints are violated
                double excess = a->radius + b->radius - dist;
                //Convert to unit vector pointing in direction of force
                vdiv_by(&displacement, vmag(&displacement));
                //Apply constants
                vmul_by(&displacement, -STERIC_FORCE_CONSTANT * excess);
                //Apply to atom a
                vadd_to(&a->force, &displacement);
            }
        }
    }
}

#define POLAR_KICK_PROB 0.0001
#define KICK_PROB 0.0003
#define WATER_RADIUS 1.4
#define DRAG_SHIELDING_DISTANCE 10.14
#define COS_DRAG_BLOCK_ANGLE 0.80901699437494745
#define KICK_VELOCITY 0.08

void water_force(struct model *m, struct steric_grid *g){
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

        struct vector displacement;
        if(do_kick){
            bool good = true;
            struct atom_list *l;
            for(l = g->cells[g->interaction_list[i]]; l; l = l->next){
                struct atom *b = m->atoms + l->atom_idx;
                if(a == b)
                    continue;

                vsub(&displacement, &kick_point, &b->position);
                if(vmag(&displacement) < b->radius + WATER_RADIUS){
                    good = false;
                    break;
                }
            }

            if(good){
                vmul_by(&kick, KICK_VELOCITY);
                vadd_to(&a->force, &kick);
            }
        }
    }
}

void drag_force(struct model *m, struct steric_grid *g){
    for(size_t i=0; i < m->num_atoms; i++){
        struct atom *a = &m->atoms[i];
        if(a->fixed)
            continue;

        struct vector displ, drag;
        bool apply = true;
        struct atom_list *l;
        for(l = g->cells[g->interaction_list[i]]; l; l = l->next){
            struct atom *b = m->atoms + l->atom_idx;
            if(a == b)
                continue;

            vsub(&displ, &b->position, &a->position);
            double dist = vmag(&displ);

            if(dist < DRAG_SHIELDING_DISTANCE){
                double dot = vdot(&displ, &a->velocity);
                if(dot / (dist * vmag(&a->velocity)) > COS_DRAG_BLOCK_ANGLE){
                    apply = false;
                    break;
                }
            }
        }

        if(apply){
            vector_copy_to(&drag, &a->velocity);
            vmul_by(&drag, m->drag_coefficient);
            vadd_to(&a->force, &drag);
        }
    }
}

