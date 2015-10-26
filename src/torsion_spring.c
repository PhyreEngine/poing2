#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "torsion_spring.h"

static void torsion_spring_force_single(struct vector *dst, struct
        torsion_spring *s, enum torsion_unit on);

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

void torsion_spring_force(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct vector *f4,
        struct torsion_spring *s){
    vector_zero(f2);
    vector_zero(f3);
    torsion_spring_force_single(f1, s, R1);
    torsion_spring_force_single(f4, s, R4);
}

void torsion_spring_force_single(struct vector *dst, struct torsion_spring *s,
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

/**
 * Here, we calculate the force on each atom of a torsion spring from a
 * potential energy function.  We use the CHARMM method (Blondel & Karplus
 * 1995: New Formulation for Derivatives of Torsion Angles and Improper
 * Torsion Angles in Molecular Mechanics: Elimination of Singularities) rather
 * than the GROMACS method (GROMACS manual, sec 4.2.13 and Allen & Tildesley
 * (pp. 330-332)).
 *
 * The force \f$\vec{F}_i\f$ on particle \f$i\f$ at coordinates \f$\vec{r_i}\f$
 * is given by
 * \f{eqnarray*}
 *      \vec{F}_i &=& -\vec_{\vec{r}}\nabla U(\phi)
 *                &=& -\frac{\partial U(\phi)}{\partial x}\vec{i}
 *                    -\frac{\partial U(\phi)}{\partial y}\vec{j}
 *                    -\frac{\partial U(\phi)}{\partial j}\vec{z}
 * \f}
 * where \f$\vec{\vec{r_i}}\f$ indicates the gradient in cartesian coordinates,
 * \f$\phi\f$ is the dihedral angle and \f$U(\phi\f$ is the potential function.
 *
 * From the chain rule,
 * \f[
 *      \frac{\partial U(\phi)}{\partial x} =
 *          \frac{\partial U(\phi)}{\partial \phi} \cdot
 *          \frac{\partial \phi}{\partial x}.
 * \f]

 * The partial derivatives \f$ \partial \phi / \partial \vec{r}_i \f$ are given
 * by Blondel and Karplus (eqns 27).
 */

void torsion_spring_force_new(
        struct vector *f1,
        struct vector *f2,
        struct vector *f3,
        struct vector *f4,
        struct torsion_spring *s){

    struct vector tmp;

    //Bond vectors according to naming in Blondel & Karplus
    struct vector F, G, H;
    vsub(&F, &s->a1->position, &s->a2->position);
    vsub(&G, &s->a2->position, &s->a3->position);
    vsub(&H, &s->a4->position, &s->a3->position);

    //Intermediate vectors, again named according to B&K
    struct vector A, B;
    vcross(&A, &F, &G);
    vcross(&B, &H, &G);

    //Get the dE/dphi part.  As a quick test, let's just use
    //-cos(phi-phi0) as the potential well. That gives us a force of
    //-sin(phi-phi0).
    //double angle = acos(vdot(&A, &B) / (vmag(&A) * vmag(&B)));
    vcross(&tmp, &B, &A);
    double angle = torsion_spring_angle(s);
    double force = -sin((angle - s->angle) / 180 * M_PI);

    //d \phi / d r_i (i.e. for first atom)
    vector_copy_to(f1, &A);
    vmul_by(f1, - force * vmag(&G) / (vmag(&A) * vmag(&A)));

    //d \phi / d r_j (i.e. for second atom)
    vector_zero(f2);
    //First term
    vector_copy_to(&tmp, &A);
    vmul_by(&tmp, force * vmag(&G) / (vmag(&A) * vmag(&A)));
    vadd_to(f2, &tmp);
    //Second term
    vector_copy_to(&tmp, &A);
    vmul_by(&tmp, force * vdot(&F, &G) / (vmag(&A) * vmag(&A) *  vmag(&G)));
    vadd_to(f2, &tmp);
    //third term
    vector_copy_to(&tmp, &B);
    vmul_by(&tmp, -force * vdot(&H, &G) / (vmag(&B) * vmag(&B) *  vmag(&G)));
    vadd_to(f2, &tmp);


    //d \phi / d r_k (i.e. for third atom)
    vector_zero(f3);
    //First term
    vector_copy_to(&tmp, &B);
    vmul_by(&tmp, force * vdot(&H, &G) / (vmag(&B) * vmag(&B) * vmag(&G)));
    vadd_to(f3, &tmp);
    //Second term
    vector_copy_to(&tmp, &A);
    vmul_by(&tmp, -force * vdot(&F, &G) / (vmag(&A) * vmag(&A) * vmag(&G)));
    vadd_to(f3, &tmp);
    //Third term
    vector_copy_to(&tmp, &B);
    vmul_by(&tmp, -force * vmag(&G) / (vmag(&B) * vmag(&B)));
    vadd_to(f3, &tmp);

    //Fourth atom
    vector_copy_to(f4, &B);
    vmul_by(f4, force * vmag(&G) / (vmag(&B) * vmag(&B)));

    //Multiply by constant
    vmul_by(f1, s->constant);
    vmul_by(f2, s->constant);
    vmul_by(f3, s->constant);
    vmul_by(f4, s->constant);
}
