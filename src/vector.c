#include <stdlib.h>
#include <math.h>
#include "vector.h"

vector vector_zero(){
    return vector_fill(0, 0, 0);
}

vector vector_copy(vector v){
    vector r = vector_zero();
    for(unsigned int i=0; i < 3; i++)
        r[i] = v[i];
    return r;
}

vector vector_fill(double x, double y, double z){
    vector v = malloc(sizeof(double) * 3);
    if(!v)
        return NULL;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
}

vector vsub_to(vector v1, vector v2){
    for(unsigned int i=0; i < 3; i++)
        v1[i] -= v2[i];
    return v1;
}

vector vadd_to(vector v1, vector v2){
    for(unsigned int i=0; i < 3; i++)
        v1[i] += v2[i];
    return v1;
}

vector vmul_by(vector v1, double s){
    for(unsigned int i=0; i < 3; i++)
        v1[i] *= s;
    return v1;
}

vector vdiv_by(vector v1, double s){
    for(unsigned int i=0; i < 3; i++)
        v1[i] /= s;
    return v1;
}

vector vsub(vector v1, vector v2){
    vector r = vector_copy(v1);
    vsub_to(r, v2);
    return r;
}

vector vadd(vector v1, vector v2){
    vector r = vector_copy(v1);
    vadd_to(r, v2);
    return r;
}

vector vmul(vector v1, double s){
    vector r = vector_copy(v1);
    vmul_by(r, s);
    return r;
}

vector vdiv(vector v1, double s){
    vector r = vector_copy(v1);
    vdiv_by(r, s);
    return r;
}

vector vcross(vector v1, vector v2){
    vector r = vector_zero();
    r[0] = v1[1]*v2[2] - v1[2]*v2[1];
    r[1] = v1[2]*v2[0] - v1[0]*v2[2];
    r[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return r;
}

double vdot(vector v1, vector v2){
    double s=0;
    for(unsigned int i=0; i < 3; i++)
        s += v1[i]*v2[i];
    return s;
}

double vmag(vector v1){
    return sqrt(vdot(v1, v1));
}
