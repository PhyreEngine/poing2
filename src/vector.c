#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "vector.h"

struct vector *vector_alloc(double x, double y, double z){
    struct vector *v = malloc(sizeof(struct vector));
    if(!v)
        return NULL;
    v->c[0] = x;
    v->c[1] = y;
    v->c[2] = z;
    return v;
}

void vector_free(struct vector *v){
    free(v);
}

void vector_zero(struct vector *v){
    for(size_t i=0; i < N; i++)
        v->c[i] = 0;
}

void vector_copy_to(struct vector *dst, struct vector *src){
    for(size_t i=0; i < N; i++)
        dst->c[i] = src->c[i];
}

void vector_fill(struct vector *v, double x, double y, double z){
    v->c[0] = x;
    v->c[1] = y;
    v->c[2] = z;
}

void vsub_to(struct vector *v1, struct vector *v2){
    vsub(v1, v1, v2);
}

void vadd_to(struct vector *v1, struct vector *v2){
    vadd(v1, v1, v2);
}

void vmul_by(struct vector *v1, double s){
    vmul(v1, v1, s);
}

void vdiv_by(struct vector *v1, double s){
    vdiv(v1, v1, s);
}

void vsub(struct vector *dst, struct vector *v1, struct vector *v2){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] - v2->c[i];
}

void vadd(struct vector *dst, struct vector *v1, struct vector *v2){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] + v2->c[i];
}

void vmul(struct vector *dst, struct vector *v1, double s){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] * s;
}

void vdiv(struct vector *dst, struct vector *v1, double s){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] / s;
}

void vcross(struct vector *dst, struct vector *v1, struct vector *v2){
    dst->c[0] = v1->c[1]*v2->c[2] - v1->c[2]*v2->c[1];
    dst->c[1] = v1->c[2]*v2->c[0] - v1->c[0]*v2->c[2];
    dst->c[2] = v1->c[0]*v2->c[1] - v1->c[1]*v2->c[0];
}

double vdot(struct vector *v1, struct vector *v2){
    double s=0;
    for(size_t i=0; i < N; i++)
        s += v1->c[i]*v2->c[i];
    return s;
}

double vmag(struct vector *v1){
    return sqrt(vdot(v1, v1));
}
