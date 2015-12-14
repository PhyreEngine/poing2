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
    s->enabled = true;
    s->cutoff = -1;

    //Default to not disabling based on handedness
    s->inner = s->outer = NULL;
    s->right_handed = false;
}

void linear_spring_free(struct linear_spring *s){
    free(s);
}

bool linear_spring_active(struct linear_spring *s){
    struct vector displacement;
    vsub(&displacement, &s->b->position, &s->a->position);
    double distance = vmag(&displacement);

    if(!s->enabled)
        return false;
    if(s->inner && s->outer){
        /*
         * Here is the situation:
         *          o--------b
         *          |
         *          |
         * a--------i
         *
         * We want to determine whether o is above or below the plane defined by
         * aib. To do this we take the cross product of ab and ai, to find the
         * normal to the plane aib. Then we take the dot product of ao and the
         * normal, and call this spring right handed if the dot product is
         * greater than zero.
         */
        struct vector ab; //Vector between A and B
        struct vector ai; //Between a and i (inner)
        struct vector oa; //Between o (outer) and a
        struct vector cross;
        vsub(&ab, &s->b->position, &s->a->position);
        vsub(&ai, &s->inner->position, &s->a->position);
        vsub(&oa, &s->outer->position, &s->a->position);
        vcross(&cross, &ab, &ai);
        double dot = vdot(&cross, &oa);
        bool rh = (dot > 0) ? true : false;
        if((rh && !s->right_handed) || (!rh && s->right_handed))
            return false;

    }
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
