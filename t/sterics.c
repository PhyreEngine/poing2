#include <stdio.h>
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

void increment(struct atom *a, struct atom *b, void *data){
    int *inc = (int*)data;
    (*inc)++;
}

int main(){
    plan(8);

    const char *springs =  "sequence = AA\n";
    struct model *m = springreader_parse_str(springs);

    vector_fill(&m->residues[0].atoms[0].position, 0, 0, 0);
    vector_fill(&m->residues[0].atoms[1].position, 1, 0, 0);
    vector_fill(&m->residues[1].atoms[0].position, 0, 0, 1);
    vector_fill(&m->residues[1].atoms[1].position, 1, 1, 1);

    struct steric_grid grid;
    steric_grid_init(&grid, 4);
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
    ok(m->residues[0].atoms[0].force.c[0] < 0, "CA1 pushed in -x");

    int num_nearby = 0;
    steric_grid_foreach_nearby(&grid, &m->residues[0].atoms[0],
            increment, &num_nearby);
    cmp_ok(num_nearby, "==", 3, "Three nearby atoms");

    done_testing();
}
