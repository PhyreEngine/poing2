#ifndef LINEAR_SPRING_H_
#define LINEAR_SPRING_H_

#include <stdbool.h>
#include "residue.h"
#include "vector.h"

#define DEFAULT_SPRING_CONSTANT 0.01

#define SC_BB_SPRING_CONSTANT 0.1

#define BB_BB_SPRING_CONSTANT 0.05

struct linear_spring {
    double distance;
    double constant;
    double cutoff;
    bool enabled;
    struct atom *a, *b;

    ///Used for determining handedness
    struct atom *inner, *outer;
    bool right_handed;
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

bool linear_spring_active(struct linear_spring *s);

void linear_spring_force(
        struct vector *dst,
        struct linear_spring *s,
        enum unit on);

#endif //LINEAR_SPRING_H_
