#include <stdlib.h>
#include "linear_spring.h"

struct linear_spring * linear_spring_alloc(double distance, double constant,
        struct atom *a, struct atom *b){
    struct linear_spring * s = malloc(sizeof(struct linear_spring));
    if(!s)
        return NULL;

    s->distance = distance;
    s->constant = constant;
    s->a = a;
    s->b = b;
    return s;
}

void linear_spring_free(struct linear_spring *s){
    free(s);
}

void linear_spring_force(
        struct vector *dst, struct linear_spring *s, enum unit on){

    struct vector displacement;
    if(on == A)
        vsub(&displacement, &s->b->position, &s->a->position);
    else
        vsub(&displacement, &s->a->position, &s->b->position);

    //Normalise displacement to get direction
    double distance = vmag(&displacement);
    vdiv_by(&displacement, distance);

    vmul(dst, &displacement, (distance - s->distance) * s->constant);
}
