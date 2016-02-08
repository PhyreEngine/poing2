#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../src/sterics.h"
#include "../src/model.h"
#include "../src/springreader.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2,
        double epsilon, const char *text){

    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

bool increment(struct atom *a, struct atom *b, void *data){
    int *inc = (int*)data;
    (*inc)++;
    return true;
}

int main(){
    plan(10);

    const char *springs =
        "[PDB]\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n"
        "ATOM      2 ALA  ALA A   1       0.000   0.000   0.000\n"
        "ATOM      3  CA  ALA A   2       0.000   0.000   0.000\n"
        "ATOM      4 ALA  ALA A   2       0.000   0.000   0.000\n"
        ;

    const char *atom_names[4] = {"N", "C", "CA", "O"};
    struct atom atoms[4];
    for(int i=0; i < 4; i++){
        atom_init(&atoms[i], i+1, atom_names[i]);
        atom_set_atom_description(&atoms[i],
                atom_description_lookup(
                    atom_names[i], strlen(atom_names[i])));
        vector_zero(&atoms[i].force);
    }
    struct model *m = model_alloc();
    m->atoms = atoms;
    m->num_atoms = 4;


    vector_fill(&atoms[0].position, 0, 0, 0);
    vector_fill(&atoms[1].position, 1, 0, 0);
    vector_fill(&atoms[2].position, 0, 0, 1);
    vector_fill(&atoms[3].position, 1, 1, 1);

    struct steric_grid grid;
    steric_grid_init(&grid, 1.5);
    steric_grid_update(&grid, m);
    is_vector(
            &grid.min,
            vector_alloc(0-GRID_BUFFER, 0-GRID_BUFFER, 0-GRID_BUFFER),
            1e-10, "Minimum grid point");
    is_vector(
            &grid.max,
            vector_alloc(1+GRID_BUFFER, 1+GRID_BUFFER, 1+GRID_BUFFER),
            1e-10, "Maximum grid point");

    steric_grid_forces(&grid, m);
    ok(atoms[0].force.c[0] < 0, "CA1 pushed in -x");

    int num_nearby = 0;
    steric_grid_foreach_nearby(&grid, &atoms[0],
            increment, &num_nearby);
    cmp_ok(num_nearby, "==", 3, "Three nearby atoms");

    //Lower cell size and see how many atoms we find
    grid.cell_size = 0.5;
    steric_grid_update(&grid, m);
    num_nearby = 0;
    steric_grid_foreach_nearby(&grid, &atoms[0],
            increment, &num_nearby);
    cmp_ok(num_nearby, "==", 0, "found 0 atoms in nearby cells");

    //Move an atom closer and ensure we find it
    atoms[1].position.c[0] = 0.2;
    num_nearby = 0;
    steric_grid_update(&grid, m);
    steric_grid_foreach_nearby(&grid, &atoms[0],
            increment, &num_nearby);
    cmp_ok(num_nearby, "==", 1, "found 1 atom in nearby cells");

    done_testing();
}
