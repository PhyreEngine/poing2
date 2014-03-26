#include <stdio.h>
#include <stdlib.h>
#include "../src/rk4.h"
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/linear_spring.h"
#include "../src/torsion_spring.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(vector v1, vector v2, const char *text){
    ok(abs(v1[0] - v2[0]) < 1e-9, "%s: x element", text);
    ok(abs(v1[1] - v2[1]) < 1e-9, "%s: y element", text);
    ok(abs(v1[2] - v2[2]) < 1e-9, "%s: z element", text);
}

int main(int argc, char **argv){
    plan(2);
/*
    struct residue *r1 = residue_alloc(AA_lookup("G", 1));
    struct residue *r2 = residue_alloc(AA_lookup("G", 1));
    struct residue residues[2];
    residues[0] = *r1;
    residues[1] = *r2;

    struct linear_spring *s = linear_spring_alloc(2.0, 1.0, r1, r2);
    vector_fill(r1->position, -0.5, 0, 0);
    vector_fill(r2->position, +0.5, 0, 0);
    struct linear_spring springs[1];
    springs[0] = *s;

    struct model *m = model_alloc();
    m->num_residues        = 2;
    m->num_linear_springs  = 1;
    m->num_torsion_springs = 0;
    m->residues            = residues;
    m->linear_springs      = springs;

    rk4_push(m, 0.1);
    ok(r1->position[0] < -0.5, "Moved r1 in the correct direction");
    ok(r2->position[0] > +0.5, "Moved r2 in the correct direction");

    for(int i=0; i < 1000; i++){
        printf("%g %g %g %g %g %g %g %g %g %g %g %g\n",
                r1->position[0], r1->position[1], r1->position[2],
                r2->position[0], r2->position[1], r2->position[2],
                r1->velocity[0], r1->velocity[1], r1->velocity[2],
                r2->velocity[0], r2->velocity[1], r2->velocity[2]);
        rk4_push(m, 0.1);
    }
    */

    struct residue *r1 = residue_alloc(AA_lookup("G", 1));
    struct residue *r2 = residue_alloc(AA_lookup("G", 1));
    struct residue *r3 = residue_alloc(AA_lookup("G", 1));
    struct residue *r4 = residue_alloc(AA_lookup("G", 1));
    struct residue residues[4];
    residues[0] = *r1;
    residues[1] = *r2;
    residues[2] = *r3;
    residues[3] = *r4;
    vector_fill(r1->position, -1, 0, 0);
    vector_fill(r2->position,  0, 0, 0);
    vector_fill(r3->position,  0, 0, 1);
    vector_fill(r4->position,  1, 0, 1);

    struct linear_spring *ls1 = linear_spring_alloc(1.0, 1.0, r1, r2);
    struct linear_spring *ls2 = linear_spring_alloc(1.0, 1.0, r2, r3);
    struct linear_spring *ls3 = linear_spring_alloc(1.0, 1.0, r3, r4);
    struct linear_spring l_springs[3];
    l_springs[0] = *ls1;
    l_springs[1] = *ls2;
    l_springs[2] = *ls3;

    struct torsion_spring *s = torsion_spring_alloc(
            r1, r2, r3, r4,
            90, .1);
    struct torsion_spring springs[1];
    springs[0] = *s;

    struct model *m = model_alloc();
    m->num_residues        = 4;
    m->num_linear_springs  = 3;
    m->num_torsion_springs = 1;
    m->residues            = residues;
    m->torsion_springs     = springs;
    m->linear_springs      = l_springs;

    for(int i=0; i < 10000; i++){
        printf("%g %g %g %g %g %g %g %g %g %g %g %g\n",
                r1->position[0], r1->position[1], r1->position[2],
                r2->position[0], r2->position[1], r2->position[2],
                r3->position[0], r3->position[1], r3->position[2],
                r4->position[0], r4->position[1], r4->position[2],
                r1->velocity[0], r1->velocity[1], r1->velocity[2],
                r2->velocity[0], r2->velocity[1], r2->velocity[2],
                r3->velocity[0], r3->velocity[1], r3->velocity[2],
                r4->velocity[0], r4->velocity[1], r4->velocity[2]);
        rk4_push(m, 0.01);
    }


    done_testing();
}
