#include <stdlib.h>
#include <math.h>
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

double torsion_spring_angle(struct torsion_spring *s){
    vector b1 = vsub(s->r2->position, s->r1->position);
    vector b2 = vsub(s->r3->position, s->r2->position);
    vector b3 = vsub(s->r4->position, s->r3->position);

    double x = vdot(
            vcross(vcross(b1, b2), vcross(b2, b3)),
            vdiv(b2, vmag(b2)));
    double y = vdot(vcross(b1, b2), vcross(b2, b3));
    return atan2(x, y) * 180 / M_PI;
}

vector torsion_spring_axis(struct torsion_spring *s){
    return vsub(s->r3->position, s->r2->position);
}

vector torsion_spring_torque(struct torsion_spring *s){
    double delta_angle = (torsion_spring_angle(s) - s->angle) / 180 * M_PI;
    vector axis = torsion_spring_axis(s);
    vdiv_by(axis, vmag(axis));
    double torque_mag = -s->constant * delta_angle;
    vmul_by(axis, torque_mag);
    return axis;
}

vector torsion_spring_force(struct torsion_spring *s, enum torsion_unit on){
    vector torque = torsion_spring_torque(s);
    vector arm = (on == R1) ?
        vsub(s->r1->position, s->r2->position)
        :
        vsub(s->r4->position, s->r3->position);
    double r = vmag(arm);
    vdiv_by(arm, r);
    return vdiv(vcross(arm, torque), (on == R1)? r : -r);
}
