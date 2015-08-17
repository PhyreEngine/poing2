#ifndef TORSION_SPRING_H_
#define TORSION_SPRING_H_

#include "residue.h"
#include "vector.h"

enum torsion_unit {R1, R4};

struct torsion_spring {
    struct atom *a1, *a2, *a3, *a4;
    double angle;
    double constant;
    double cutoff;
    bool enabled;
};

struct torsion_spring * torsion_spring_alloc(
        struct atom *a1, struct atom *a2,
        struct atom *a3, struct atom *a4,
        double angle, double constant);

void torsion_spring_init(
        struct torsion_spring *s,
        struct atom *a1, struct atom *a2,
        struct atom *a3, struct atom *a4,
        double angle, double constant);
void torsion_spring_free(struct torsion_spring *s);
void torsion_spring_axis(struct vector *dst, struct torsion_spring *s);
void torsion_spring_torque(struct vector *dst, struct torsion_spring *s);
double torsion_spring_angle(struct torsion_spring *s);

void torsion_spring_force(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct vector *f4,
        struct torsion_spring *s);
void torsion_spring_force_new(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct vector *f4,
        struct torsion_spring *s);

#endif /* TORSION_SPRING_H_ */

