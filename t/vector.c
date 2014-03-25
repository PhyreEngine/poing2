#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/vector.h"
#include "tap.h"

void is_vector(vector v1, vector v2, const char *text){
    ok(abs(v1[0] - v2[0]) < 1e-9, "%s: x element", text);
    ok(abs(v1[1] - v2[1]) < 1e-9, "%s: y element", text);
    ok(abs(v1[2] - v2[2]) < 1e-9, "%s: z element", text);
}

int main(int argc, char **argv){
    plan(32);

    vector f = vector_zero();
    vector_fill(f, 3, 2, 1);
    is_vector(f, vector_create(3, 2, 1), "Filled vector");

    vector v1 = vector_create(1, 2, 3);
    vector v2 = vector_create(3, 2, 1);

    ok(abs(vmag(v1) - sqrt(14)) < 1e-9, "Vector magnitude");

    vector diff  = vsub(v1, v2);
    vector sum   = vadd(v1, v2);
    vector cross = vcross(v1, v2);
    double dot   = vdot(v1, v2);
    vector mul   = vmul(v1, 3.);
    vector div   = vdiv(v1, 3.);

    ok(abs(dot - 10) < 1e-9, "Dot product");
    is_vector(diff,  vector_create(-2,    0,     2),   "Subtraction");
    is_vector(sum,   vector_create(4,     4,     4),   "Addition");
    is_vector(cross, vector_create(-4,    8,     -4),  "Cross product");
    is_vector(mul,   vector_create(3,     6,     9),   "Scalar multiplication");
    is_vector(div,   vector_create(1.0/3, 2.0/3, 1.0), "Scalar division");

    vsub_to(v1, v2);
    is_vector(v1,  vector_create(-2, 0, 2),   "Subtraction side effect");
    vadd_to(v1, v2);
    is_vector(v1,  vector_create(1, 2, 3),   "Addition side effect");
    vmul_by(v1, 3.);
    is_vector(v1,  vector_create(3, 6, 9),   "Multiplication side effect");
    vdiv_by(v1, 3.);
    is_vector(v1,  vector_create(1, 2, 3),   "Division side effect");

    done_testing();
}
