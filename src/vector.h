#ifndef VECTOR_H_
#define VECTOR_H_

typedef double* vector;

vector vector_zero();
vector vector_copy(vector v);
vector vector_fill(double x, double y, double z);
vector vsub_to(vector v1, vector v2);
vector vadd_to(vector v1, vector v2);
vector vmul_by(vector v1, double s);
vector vdiv_by(vector v1, double s);
vector vsub(vector v1, vector v2);
vector vadd(vector v1, vector v2);
vector vmul(vector v1, double s);
vector vdiv(vector v1, double s);
vector vcross(vector v1, vector v2);
double vdot(vector v1, vector v2);
double vmag(vector v1);

#endif /* VECTOR_H_ */

