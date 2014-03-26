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

    //Bond vectors
    struct vector b1, b2, b3;
    vsub(&b1, &s->r2->position, &s->r1->position);
    vsub(&b2, &s->r3->position, &s->r2->position);
    vsub(&b3, &s->r4->position, &s->r3->position);

    struct vector cross_b1_b2, cross_b2_b3;
    vcross(&cross_b1_b2, &b1, &b2);
    vcross(&cross_b2_b3, &b2, &b3);

    struct vector cross1;
    vcross(&cross1, &cross_b1_b2, &cross_b2_b3);

    struct vector norm_b2;
    double b2_mag = vmag(&b2);
    vdiv(&norm_b2, &b2, b2_mag);

    double y = vdot(&cross1, &norm_b2);
    double x = vdot(&cross_b1_b2, &cross_b2_b3);
    return atan2(y, x) * 180 / M_PI;
}

void torsion_spring_axis(struct vector *dst, struct torsion_spring *s){
    vsub(dst, &s->r3->position, &s->r2->position);
}

void torsion_spring_torque(struct vector *dst, struct torsion_spring *s){
    double delta_angle = (torsion_spring_angle(s) - s->angle) / 180 * M_PI;

    //Normalised axis
    struct vector axis;
    torsion_spring_axis(&axis, s);
    vdiv_by(&axis, vmag(&axis));

    //Magnitude of torque
    double torque_mag = -s->constant * delta_angle;

    //Combine direction and magnitude
    vmul(dst, &axis, torque_mag);
}

void torsion_spring_force(struct vector *dst, struct torsion_spring *s,
        enum torsion_unit on){

    struct vector torque, arm;
    torsion_spring_torque(&torque, s);
    if(on == R1)
        vsub(&arm, &s->r1->position, &s->r2->position);
    else
        vsub(&arm, &s->r4->position, &s->r3->position);

    //Direction of arm
    double r   = vmag(&arm);
    double tau = vmag(&torque);
    vdiv_by(&arm, r);
    vdiv_by(&torque, tau);

    //Direction of force
    vcross(dst, &torque, &arm);
    //Correct magnitude by arm length
    vmul_by(dst, (on == R1) ? -tau/r : tau/r);
}
