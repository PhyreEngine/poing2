#include <stdlib.h>
#include <stdbool.h>
#include "vector.h"

struct vector *vector_alloc(double x, double y, double z){
    struct vector *v = malloc(sizeof(struct vector));
    if(!v)
        return NULL;
    v->c[0] = x;
    v->c[1] = y;
    v->c[2] = z;
    return v;
}

void vector_free(struct vector *v){
    free(v);
}

extern inline void vector_zero(struct vector *);
extern inline void vector_copy_to(struct vector *dst, struct vector *src);
extern inline void vector_fill(struct vector *v, double x, double y, double z);
extern inline void vsub_to(struct vector *v1, struct vector *v2);
extern inline void vadd_to(struct vector *v1, struct vector *v2);
extern inline void vmul_by(struct vector *v1, double s);
extern inline void vdiv_by(struct vector *v1, double s);
extern inline void vsub(struct vector *dst, struct vector *v1, struct vector *v2);
extern inline void vadd(struct vector *dst, struct vector *v1, struct vector *v2);
extern inline void vmul(struct vector *dst, struct vector *v1, double s);
extern inline void vdiv(struct vector *dst, struct vector *v1, double s);
extern inline void vcross(struct vector *dst, struct vector *v1, struct vector *v2);
extern inline double vdot(struct vector *v1, struct vector *v2);
extern inline double vmag(struct vector *v1);
extern inline double vmag_sq(struct vector *v1);
extern inline void vector_rand(struct vector *dst, double min_phi, double max_phi);
extern inline void vector_spherical_coords(struct vector *dst, struct vector *v);
extern inline void vrot_x(struct vector *dst, struct vector *v, double theta);
extern inline void vrot_y(struct vector *dst, struct vector *v, double theta);
extern inline void vrot_z(struct vector *dst, struct vector *v, double theta);
extern inline void vrot_axis(struct vector *dst, struct vector *axis, struct vector *v,
        double theta);

extern inline void vmin_elems(struct vector *dst, struct vector *v);
extern inline void vmax_elems(struct vector *dst, struct vector *v);
