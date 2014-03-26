#ifndef VECTOR_H_
#define VECTOR_H_

#define N 3
typedef double real;

struct vector { real c[N]; };

struct vector * vector_alloc(double x, double y, double z);
void vector_free(struct vector *v);

void vector_zero(struct vector *);
void vector_copy_to(struct vector *dst, struct vector *src);
void vector_fill(struct vector *v, double x, double y, double z);
void vsub_to(struct vector *v1, struct vector *v2);
void vadd_to(struct vector *v1, struct vector *v2);
void vmul_by(struct vector *v1, double s);
void vdiv_by(struct vector *v1, double s);
void vsub(struct vector *dst, struct vector *v1, struct vector *v2);
void vadd(struct vector *dst, struct vector *v1, struct vector *v2);
void vmul(struct vector *dst, struct vector *v1, double s);
void vdiv(struct vector *dst, struct vector *v1, double s);
void vcross(struct vector *dst, struct vector *v1, struct vector *v2);
double vdot(struct vector *v1, struct vector *v2);
double vmag(struct vector *v1);

#endif /* VECTOR_H_ */

