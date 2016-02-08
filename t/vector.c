#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2, double epsilon, const char *text){
    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(){
    plan(62);

    struct vector zero = {.c = {0.0, 0.0, 0.0}};
    struct vector cmp  = {.c = {1.0, 2.0, 3.0}};
    struct vector v1, v2, v3, result;

    vector_zero(&v1);
    is_vector(&v1, &zero, 1e-10, "Zeroed vector");

    vector_copy_to(&v1, &cmp);
    is_vector(&v1, &cmp, 1e-10, "Copied vector");

    vector_fill(&v2, 1.0, 2.0, 3.0);
    is_vector(&v2, &cmp, 1e-10, "Filled vector");

    vector_fill(&v1, 0.0, 0.0, 0.0);
    vector_fill(&result, -1.0, -2.0, -3.0);
    vsub_to(&v1, &cmp);
    is_vector(&v1, &result, 1e-10, "In-place subtraction");

    vector_fill(&v1, 2.0, 4.0, 6.0);
    vector_fill(&result, 1.0, 2.0, 3.0);
    vsub(&v2, &v1, &cmp);
    is_vector(&v2, &result, 1e-10, "Subtraction");

    vector_fill(&v1, 1.0, 2.0, 3.0);
    vector_fill(&result, 2.0, 4.0, 6.0);
    vadd_to(&v1, &cmp);
    is_vector(&v1, &result, 1e-10, "In-place addition");

    vector_fill(&v1, 1.0, 2.0, 3.0);
    vector_fill(&result, 2.0, 4.0, 6.0);
    vadd(&v2, &v1, &cmp);
    is_vector(&v2, &result, 1e-10, "Addition");

    vector_fill(&v1, 1.0, 2.0, 3.0);
    vector_fill(&result, 2.0, 4.0, 6.0);
    vmul_by(&v1, 2);
    is_vector(&v1, &result, 1e-10, "In-place scalar multiplication");

    vector_fill(&v1, 1.0, 2.0, 3.0);
    vector_fill(&result, 2.0, 4.0, 6.0);
    vmul(&v2, &v1, 2);
    is_vector(&v2, &result, 1e-10, "Scalar multiplication");

    vector_fill(&v1, 2.0, 4.0, 6.0);
    vector_fill(&result, 1.0, 2.0, 3.0);
    vdiv_by(&v1, 2);
    is_vector(&v1, &result, 1e-10, "In-place scalar division");

    vector_fill(&v1, 2.0, 4.0, 6.0);
    vector_fill(&result, 1.0, 2.0, 3.0);
    vdiv(&v2, &v1, 2);
    is_vector(&v2, &result, 1e-10, "Scalar division");

    vector_fill(&v1, 1.0, 2.0, 3.0);
    vector_fill(&v2, 3.0, 2.0, 1.0);
    vector_fill(&result, -4.0, 8.0, -4.0);

    vcross(&v3, &v1, &v2);
    is_vector(&v3, &result, 1e-10, "Cross product");

    double dot = vdot(&v1, &v2);
    fis(dot, 10.0, 1e-10, "Dot product");

    double mag = vmag(&v1);
    fis(mag*mag, 14.0, 1e-10, "Vector magnitude");

    vector_fill(&result, 1.0, 2.0, 3.0);
    struct vector *alloced = vector_alloc(1.0, 2.0, 3.0);
    is_vector(alloced, &result, 1e-10, "Allocated vector");
    vector_free(alloced);

    //Rotation of vectors
    vector_fill(&v1, 0, 0, 1);
    vector_fill(&result, 0, -1, 0);
    vrot_x(&v2, &v1, M_PI/2);
    is_vector(&v2, &result, 1e-10, "Rotation about x axis");

    vector_fill(&v1, 0, 0, 1);
    vector_fill(&result, 1, 0, 0);
    vrot_y(&v2, &v1, M_PI/2);
    is_vector(&v2, &result, 1e-10, "Rotation about y axis");

    vector_fill(&v1, 1, 0, 0);
    vector_fill(&result, 0, 1, 0);
    vrot_z(&v2, &v1, M_PI/2);
    is_vector(&v2, &result, 1e-10, "Rotation about z axis");

    vector_fill(&v1, 1, 0, 0);
    vector_fill(&v2, 0, 0, 1);
    vector_fill(&result, 0, 1, 0);
    vrot_axis(&v3, &v2, &v1, M_PI/2);
    is_vector(&v3, &result, 1e-10, "Rotation about arbitrary axis");

    //Spherical coords
    vector_fill(&v1, 0, 0, 1);
    vector_fill(&result, 1, 0, 0);
    vector_spherical_coords(&v2, &v1);
    is_vector(&v2, &result, 1e-10, "Spherical coordinates (0, 0, 1)");

    vector_fill(&v1, 1, 0, 0);
    vector_fill(&result, 1, 0, M_PI/2);
    vector_spherical_coords(&v2, &v1);
    is_vector(&v2, &result, 1e-10, "Spherical coordinates (1, 0, 0)");

    vector_fill(&v1, 1, 1, 1);
    vector_fill(&result, sqrt(3), M_PI/4, acos(1/sqrt(3)));
    vector_spherical_coords(&v2, &v1);
    is_vector(&v2, &result, 1e-10, "Spherical coordinates (1, 1, 1)");

    done_testing();
}
