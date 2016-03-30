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
    plan(9);

    struct atom atoms[4];
    for(int i=0; i < 4; i++){
        vector_zero(&atoms[i].force);
        atoms[i].radius = 1;
    }
    struct model *m = model_alloc();
    m->atoms = atoms;
    m->num_atoms = 4;

    vector_fill(&atoms[0].position, -1, 0, 0);
    vector_fill(&atoms[1].position, 1, -1, 0);
    vector_fill(&atoms[2].position, 0, 0, 1);
    vector_fill(&atoms[3].position, 1, 1, -1);

    //Check that the origin makes sense with these atoms
    struct steric_grid grid;
    steric_grid_init(&grid, 10, 5, m->num_atoms);
    steric_grid_find_origin(&grid, m);
    ok(grid.origin.c[0] < -1.0, "Grid origin x");
    ok(grid.origin.c[1] < -1.0, "Grid origin y");
    ok(grid.origin.c[2] < -1.0, "Grid origin z");

    //Now try some basic force tests
    m->num_atoms = 2;
    vector_fill(&atoms[0].position, 0, 0, 0);
    vector_fill(&atoms[1].position, 1, 0, 0);

    steric_grid_update(&grid, m);
    steric_grid_build_ilists(&grid, m);
    steric_grid_forces(&grid, m);

    ok(atoms[0].force.c[0] < 0, "Forced into the -x direction");
    ok(atoms[1].force.c[0] > 0, "Forced into the +x direction");
    fis(atoms[0].force.c[1], 0, 1e-6, "Y component = 0");
    fis(atoms[0].force.c[2], 0, 1e-6, "Z component = 0");
    fis(atoms[1].force.c[1], 0, 1e-6, "Y component = 0");
    fis(atoms[1].force.c[2], 0, 1e-6, "Z component = 0");

    done_testing();
}
