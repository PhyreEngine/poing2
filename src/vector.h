#ifndef VECTOR_H_
#define VECTOR_H_

#include <math.h>

#define N 3
typedef double real;

struct vector { real c[N]; };

struct vector * vector_alloc(double x, double y, double z);
void vector_free(struct vector *v);

inline void vector_zero(struct vector *v){
    for(size_t i=0; i < N; i++)
        v->c[i] = 0;
}

inline void vector_copy_to(struct vector *dst, struct vector *src){
    for(size_t i=0; i < N; i++)
        dst->c[i] = src->c[i];
}


inline void vector_fill(struct vector *v, double x, double y, double z){
    v->c[0] = x;
    v->c[1] = y;
    v->c[2] = z;
}

inline void vsub(struct vector *dst, struct vector *v1, struct vector *v2){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] - v2->c[i];
}

inline void vadd(struct vector *dst, struct vector *v1, struct vector *v2){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] + v2->c[i];
}

inline void vmul(struct vector *dst, struct vector *v1, double s){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] * s;
}


inline void vdiv(struct vector *dst, struct vector *v1, double s){
    for(size_t i=0; i < N; i++)
        dst->c[i] = v1->c[i] / s;
}

inline void vsub_to(struct vector *v1, struct vector *v2){
    vsub(v1, v1, v2);
}

inline void vadd_to(struct vector *v1, struct vector *v2){
    vadd(v1, v1, v2);
}

inline void vmul_by(struct vector *v1, double s){
    vmul(v1, v1, s);
}

inline void vdiv_by(struct vector *v1, double s){
    vdiv(v1, v1, s);
}


inline void vcross(struct vector *dst, struct vector *v1, struct vector *v2){
    dst->c[0] = v1->c[1]*v2->c[2] - v1->c[2]*v2->c[1];
    dst->c[1] = v1->c[2]*v2->c[0] - v1->c[0]*v2->c[2];
    dst->c[2] = v1->c[0]*v2->c[1] - v1->c[1]*v2->c[0];
}

inline double vdot(struct vector *v1, struct vector *v2){
    double s=0;
    for(size_t i=0; i < N; i++)
        s += v1->c[i]*v2->c[i];
    return s;
}

inline double vmag(struct vector *v1){
    return sqrt(vdot(v1, v1));
}

/**
 * Find a random vector within a given azimuthal angle.
 *
 * This produces a random vector with an azimuthal angle between min_phi and
 * max_phi.
 */

inline void vector_rand(struct vector *dst, double min_phi, double max_phi){
    double cos_max_phi = cos(max_phi);
    double cos_min_phi = cos(min_phi);

    double i = (double)rand() / RAND_MAX;
    double j = (double)rand() / RAND_MAX;
    double z   = cos_min_phi - i * (cos_min_phi - cos_max_phi);
    double phi = j * M_PI * 2;

    dst->c[0] = sqrt(1 - z*z) * cos(phi);
    dst->c[1] = sqrt(1 - z*z) * sin(phi);
    dst->c[2] = z;
}

inline void vector_spherical_coords(struct vector *dst, struct vector *v){
    dst->c[0] = vmag(v);
    dst->c[1] = atan2(v->c[1], v->c[0]);
    dst->c[2] = acos(v->c[2] / dst->c[0]);
}

inline void vrot_x(struct vector *dst, struct vector *v, double theta){
    dst->c[0] = v->c[0];
    dst->c[1] = cos(theta) * v->c[1] - sin(theta) * v->c[2];
    dst->c[2] = sin(theta) * v->c[1] + cos(theta) * v->c[2];
}

inline void vrot_y(struct vector *dst, struct vector *v, double theta){
    dst->c[0] =  cos(theta) * v->c[0] + sin(theta) * v->c[2];
    dst->c[1] = v->c[1];
    dst->c[2] = -sin(theta) * v->c[0] + cos(theta) * v->c[2];
}

inline void vrot_z(struct vector *dst, struct vector *v, double theta){
    dst->c[0] = cos(theta) * v->c[0] - sin(theta) * v->c[1];
    dst->c[1] = sin(theta) * v->c[0] + cos(theta) * v->c[1];
    dst->c[2] = v->c[2];
}

inline void vrot_axis(struct vector *dst, struct vector *axis, struct vector *v,
        double theta){
    double ux = axis->c[0];
    double uy = axis->c[1];
    double uz = axis->c[2];
    double ct = cos(theta);
    double st = sin(theta);

    dst->c[0] = v->c[0] * (ct + ux*ux*(1-ct))
        + v->c[1] * (ux*uy*(1-ct) - uz*st)
        + v->c[2] * (ux*uz*(1-ct) + uy*st);

    dst->c[1] = v->c[0] * (uy*ux*(1-ct) + uz*st)
        + v->c[1] * (ct + uy*uy*(1-ct))
        + v->c[2] * (uy*uz*(1-ct) - ux*st);

    dst->c[2] = v->c[0] * (uz*ux*(1-ct) - uy*st)
        + v->c[1] * (uz*uy*(1-ct) + ux*st)
        + v->c[2] * (ct + uz*uz*(1-ct));
}

inline void vmin_elems(struct vector *dst, struct vector *v){
    for(size_t i=0; i < N; i++)
        if(v->c[i] < dst->c[i])
            dst->c[i] = v->c[i];
}

inline void vmax_elems(struct vector *dst, struct vector *v){
    for(size_t i=0; i < N; i++)
        if(v->c[i] > dst->c[i])
            dst->c[i] = v->c[i];
}

#endif /* VECTOR_H_ */

