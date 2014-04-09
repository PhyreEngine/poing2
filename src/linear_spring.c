#include <stdlib.h>
#include <math.h>
#include "linear_spring.h"

struct linear_spring * linear_spring_alloc(double distance, double constant,
        struct atom *a, struct atom *b){
    struct linear_spring * s = malloc(sizeof(struct linear_spring));
    if(!s)
        return NULL;

    linear_spring_init(s, distance, constant, a, b);
    return s;
}

void linear_spring_init(struct linear_spring *s,
        double distance, double constant,
        struct atom *a, struct atom *b){
    s->distance = distance;
    s->constant = constant;
    s->a = a;
    s->b = b;
    s->cutoff = -1;
}

void linear_spring_free(struct linear_spring *s){
    free(s);
}

bool linear_spring_active(struct linear_spring *s){
    struct vector displacement;
    vsub(&displacement, &s->b->position, &s->a->position);
    double distance = vmag(&displacement);

    if(s->cutoff < 0 || fabs(distance - s->distance) < s->cutoff)
        return true;
    return false;
}

void linear_spring_force(
        struct vector *dst, struct linear_spring *s, enum unit on){

    struct vector displacement;
    //Don't apply if outside the cutoffs
    if(!linear_spring_active(s)){
        vector_zero(dst);
        return;
    }

    if(on == A)
        vsub(&displacement, &s->b->position, &s->a->position);
    else
        vsub(&displacement, &s->a->position, &s->b->position);

    //Normalise displacement to get direction
    double distance = vmag(&displacement);
    vdiv_by(&displacement, distance);

    vmul(dst, &displacement, (distance - s->distance) * s->constant);
}
