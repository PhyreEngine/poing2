#include <stdlib.h>
#include "linear_spring.h"

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
