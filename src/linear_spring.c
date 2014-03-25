#include <stdlib.h>
#include "linear_spring.h"
#include "vector.h"

struct linear_spring * linear_spring_alloc(double distance, double constant,
        struct residue *a, struct residue *b){
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

vector linear_spring_force(struct linear_spring *s, enum unit on){
    vector displacement = (on == A) ?
        vsub(s->b->position, s->a->position)
        :
        vsub(s->a->position, s->b->position);

    double distance = vmag(displacement);
    vector direction = vdiv(displacement, distance);

    vector force = vmul(direction, (distance - s->distance) * s->constant);
    return force;
}
