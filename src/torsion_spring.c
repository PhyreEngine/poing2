#include <stdlib.h>
#include <math.h>
#include "torsion_spring.h"

struct torsion_spring * torsion_spring_alloc(
        struct atom *a1, struct atom *a2,
        struct atom *a3, struct atom *a4,
        double angle, double constant){

    struct torsion_spring *s = malloc(sizeof(struct torsion_spring));
    if(!s)
        return NULL;
    torsion_spring_init(s, a1, a2, a3, a4, angle, constant);
    return s;
}

void torsion_spring_init(struct torsion_spring *s,
        struct atom *a1, struct atom *a2,
        struct atom *a3, struct atom *a4,
        double angle, double constant){
    s->a1 = a1;
    s->a2 = a2;
    s->a3 = a3;
    s->a4 = a4;
    s->enabled = true;
    s->angle = angle;
    s->constant = constant;
    s->cutoff = -1;
}

void torsion_spring_free(struct torsion_spring *s){
    free(s);
}

double torsion_spring_angle(struct torsion_spring *s){

    //Bond vectors
    struct vector b1, b2, b3;
    vsub(&b1, &s->a2->position, &s->a1->position);
    vsub(&b2, &s->a3->position, &s->a2->position);
    vsub(&b3, &s->a4->position, &s->a3->position);

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
    vsub(dst, &s->a3->position, &s->a2->position);
}

void torsion_spring_torque(struct vector *dst, struct torsion_spring *s){
    double delta_angle = (torsion_spring_angle(s) - s->angle);
    if(delta_angle < -180)
        delta_angle += 360;
    else if(delta_angle > 180)
        delta_angle -= 360;
    delta_angle *= M_PI / 180;

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

    double angle = torsion_spring_angle(s);
    if(!s->enabled || (s->cutoff > 0 && fabs(angle - s->angle) > s->cutoff)){
        vector_zero(dst);
        return;
    }

    struct vector torque, arm;
    torsion_spring_torque(&torque, s);
    if(on == R1)
        vsub(&arm, &s->a1->position, &s->a2->position);
    else
        vsub(&arm, &s->a4->position, &s->a3->position);

    /*
     * We want to project the arm to a vector perpendicular to the axis,
     * because the force should only depend on this perpendicular distance to
     * the axis.
     *
     *     r4 /
     *       /
     *      /
     *  r3 /___ <- we want this vector (a)
     *     |
     *     |
     *   __|r2
     *    /
     *   /r1
     *
     * We can get vector a by first finding the projection of (r4-r3) onto the
     * axis, p(r4-r3). Then a is (r4-r3) - p(r4-r3).
     *
     */

    struct vector axis;
    torsion_spring_axis(&axis, s);
    double axis_mag = vmag(&axis);
    double proj_mag = vdot(&arm, &axis);

    struct vector proj;
    vmul(&proj, &axis, proj_mag / axis_mag / axis_mag);
    struct vector rej;
    vsub(&rej, &arm, &proj);

    //Direction of arm
    double r   = vmag(&rej);
    double tau = vmag(&torque);
    if(r == 0 || tau == 0)
        return;
    vdiv_by(&rej, r);
    vdiv_by(&torque, tau);

    //Direction of force
    vcross(dst, &torque, &rej);
    //Correct magnitude by arm length
    vmul_by(dst, (on == R1) ? -tau/r : tau/r);
}
