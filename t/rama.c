#include <stdio.h>
#include <stdlib.h>
//#include "../src/rama.h"
#include "../src/vector.h"
#include "tap.h"

void is_vector(struct vector *v1, struct vector *v2,
        double epsilon, const char *text){

    for(size_t i=0; i < N; i++)
        fis(v1->c[i], v2->c[i], epsilon, "%s: element %d", text, i);
}

int main(){
    plan(30);

    const char *springs =
        "[PDB]\n"
        "ATOM      1  O   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      2  N   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      3  C   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      4  CA  ALA A   1       0.000   0.000   0.000\n"
        "ATOM      5  ALA ALA A   1       0.000   0.000   0.000\n"
        "ATOM      5  O   ALA A   2       0.000   0.000   0.000\n"
        "ATOM      6  N   ALA A   2       0.000   0.000   0.000\n"
        "ATOM      7  C   ALA A   2       0.000   0.000   0.000\n"
        "ATOM      8  CA  ALA A   2       0.000   0.000   0.000\n"
        "ATOM      9  ALA ALA A   2       0.000   0.000   0.000\n"
        ;
//    struct model *m = springreader_parse_str(springs);
//    struct rama *r  = rama_init(&m->residues[0], &m->residues[1]);

    done_testing();
}
