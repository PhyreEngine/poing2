#ifndef TORSION_SPRING_H_
#define TORSION_SPRING_H_

#include "residue.h"
#include "vector.h"

struct torsion_spring {
    struct residue *r1, *r2, *r3, *r4;
    double angle;
    double constant;
};

struct torsion_spring * torsion_spring_alloc(
        struct residue *r1, struct residue *r2,
        struct residue *r3, struct residue *r4,
        double angle, double constant);
void torsion_spring_free(struct torsion_spring *s);
vector torsion_spring_axis(struct torsion_spring *s);
vector torsion_spring_torque(struct torsion_spring *s);
double torsion_spring_angle(struct torsion_spring *s);
vector torsion_spring_force(struct torsion_spring *s);

#endif /* TORSION_SPRING_H_ */

