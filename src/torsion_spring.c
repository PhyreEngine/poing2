#include <stdlib.h>
#include "torsion_spring.h"

struct torsion_spring * torsion_spring_alloc(
        struct residue *r1, struct residue *r2,
        struct residue *r3, struct residue *r4,
        double angle, double constant){

    struct torsion_spring *s = malloc(sizeof(struct torsion_spring));
    if(!s)
        return NULL;
    s->r1 = r1;
    s->r2 = r2;
    s->r3 = r3;
    s->r4 = r4;
    s->angle = angle;
    s->constant = constant;
    return s;
}
void torsion_spring_free(struct torsion_spring *s){
    free(s);
}
