#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>

struct model;
struct linear_spring;
struct torsion_spring;
struct bond_angle_spring;

struct model_debug {
    ///Files to print to
    FILE *linear, *angle, *torsion, *rama;
    double interval;
};

void debug_begin(struct model *m);
void debug_linear(struct model *m, struct linear_spring *s);
void debug_torsion(struct model *m, struct torsion_spring *s);
void debug_angle(struct model *m, struct bond_angle_spring *s);

#endif /* DEBUG_H */
