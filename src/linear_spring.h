#ifndef LINEAR_SPRING_H_
#define LINEAR_SPRING_H_

#include "residue.h"
#include "vector.h"

struct linear_spring {
    double distance;
    double constant;
    struct atom *a, *b;
};

enum unit { A, B };

struct linear_spring * linear_spring_alloc(
        double distance, double constant,
        struct atom *a, struct atom *b);

void linear_spring_init(
        struct linear_spring *s,
        double distance, double constant,
        struct atom *a, struct atom *b);

void linear_spring_free(struct linear_spring *s);

void linear_spring_force(
        struct vector *dst,
        struct linear_spring *s,
        enum unit on);

#endif //LINEAR_SPRING_H_
