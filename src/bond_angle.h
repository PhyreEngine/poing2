#ifndef BOND_ANGLE_H_
#define BOND_ANGLE_H_

#include <stdbool.h>
#include "residue.h"
#include "vector.h"

struct bond_angle_spring {
    struct atom *a1, *a2, *a3;
    double angle;
    double constant;
    double cutoff;
    bool enabled;
};

struct bond_angle_spring * bond_angle_spring_alloc(
        struct atom *a1, struct atom *a2, struct atom *a3,
        double angle, double constant);

void bond_angle_spring_init(
        struct bond_angle_spring *s,
        struct atom *a1, struct atom *a2, struct atom *a3,
        double angle, double constant);
void bond_angle_spring_free(struct bond_angle_spring *s);

void bond_angle_force(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct bond_angle_spring *s);

#endif /* BOND_ANGLE_H_ */
