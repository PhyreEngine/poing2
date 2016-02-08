#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/model.h"
#include "../src/record.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2,
        double epsilon, const char *text){
    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(){
    plan(6);

    const int natoms = 3;
    struct atom atoms[natoms];
    struct model m = {.num_atoms = 3, .atoms = atoms};

    for(size_t i=0; i < natoms; i++){
        vector_zero(&atoms[i].position);
        atoms[i].fixed = false;
    }

    struct record pos;
    record_init(&pos, &m, 11);

    for(size_t i=0; i < 10; i++){
        for(size_t j=0; j < natoms; j++)
            atoms[j].position.c[0] += 0.1;
        record_add(&pos, &m);
    }

    fis(pos.avg_jitter[0], 0.1, 1e-6, "Avg = 0.1");
    fis(pos.avg_jitter[1], 0.1, 1e-6, "Avg = 0.1");
    fis(pos.avg_jitter[2], 0.1, 1e-6, "Avg = 0.1");

    for(size_t i=0; i < 5; i++){
        for(size_t j=0; j < natoms; j++)
            atoms[j].position.c[0] += 0.2;
        record_add(&pos, &m);
    }

    fis(pos.avg_jitter[0], 0.15, 1e-6, "Avg = 0.15");
    fis(pos.avg_jitter[1], 0.15, 1e-6, "Avg = 0.15");
    fis(pos.avg_jitter[2], 0.15, 1e-6, "Avg = 0.15");

    done_testing();
}
