#ifndef LINEAR_SPRING_H_
#define LINEAR_SPRING_H_

#include "residue.h"

struct linear_spring {
    double distance;
    double constant;
    struct residue *a, *b;
};

enum unit { A, B };

struct linear_spring * linear_spring_alloc(double distance, double constant,
        struct residue *a, struct residue *b);
void linear_spring_free(struct linear_spring *s);

vector linear_spring_force(struct linear_spring *s, enum unit on);

#endif //LINEAR_SPRING_H_
