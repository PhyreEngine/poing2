#include <stdio.h>
#include <stdlib.h>
#include "../src/linear_spring.h"
#include "../src/residue.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(vector v1, vector v2, const char *text){
    ok(abs(v1[0] - v2[0]) < 1e-9, "%s: x element", text);
    ok(abs(v1[1] - v2[1]) < 1e-9, "%s: y element", text);
    ok(abs(v1[2] - v2[2]) < 1e-9, "%s: z element", text);
}

int main(int argc, char **argv){
    plan(24);
    struct residue *a, *b;
    struct linear_spring *s;

    a = residue_alloc(AA_lookup("G", 1));
    b = residue_alloc(AA_lookup("G", 1));

    a->position = vector_fill(-0.5, 0, 0);
    b->position = vector_fill( 0.5, 0, 0);
    s = linear_spring_alloc(1, 1.0, a, b);

    is_vector(linear_spring_force(s, A), vector_zero(), "Force = 0 ");
    is_vector(linear_spring_force(s, B), vector_zero(), "Force = 0 ");

    s->distance = 2.0;
    is_vector(linear_spring_force(s, A), vector_fill(-1, 0, 0), "X: -1 ");
    is_vector(linear_spring_force(s, B), vector_fill( 1, 0, 0), "X:  1 ");

    s->distance = 0.5;
    is_vector(linear_spring_force(s, A), vector_fill( .5, 0, 0), "X:  0.5");
    is_vector(linear_spring_force(s, B), vector_fill(-.5, 0, 0), "X: -0.5 ");

    s->distance = 2.0;
    s->constant = 2;
    is_vector(linear_spring_force(s, A), vector_fill(-2, 0, 0), "X: -2 ");
    is_vector(linear_spring_force(s, B), vector_fill( 2, 0, 0), "X:  2 ");
    done_testing();
}
