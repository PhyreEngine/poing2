#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/springreader.h"
#include "../src/model.h"
#include "tap.h"

const char *pdb =
"timestep = 0.01\n"
"synth_time = 100\n"
"drag_coefficient = 0.1\n"
"#List of springs (i, j, distance)\n"
"[PDB]\n"
"ATOM      1  CA  GLU A   1       0.000   0.000   0.000\n"
"ATOM      2 GLU  GLU A   1       0.000   0.000   0.000\n"
"ATOM      3  CA  VAL A   2       0.000   0.000   0.000\n"
"ATOM      4 VAL  VAL A   2       0.000   0.000   0.000\n"
"ATOM      5  CA  TYR A   3       0.000   0.000   0.000\n"
"ATOM      6 TYR  TYR A   3       0.000   0.000   0.000\n"
"ATOM      7  CA  LEU A   4       0.000   0.000   0.000\n"
"ATOM      8 TYR  LEU A   4       0.000   0.000   0.000\n"
"[Linear]\n"
"1 CA 2 CA 4.0\n"
"1 CA 3 CA 3.0 0.2\n"
"3 CA 4 CA 4.0 0.3 0.1\n"
"[Torsion]\n"
"1 CA 2 CA 3 CA 4 TYR 40 0.1\n";

void check_model(struct model *m){
    fis(m->timestep, 0.01, 1e-10, "Timestep");
    fis(m->synth_time, 100, 1e-10, "Synth time");
    fis(m->drag_coefficient, 0.1, 1e-10, "Drag coefficient");

    //3 springs
    cmp_ok(m->num_linear_springs, "==",  3, "Read three linear springs");
    cmp_ok(m->num_torsion_springs, "==", 1, "Read one torsion spring");
    cmp_ok(m->num_residues, "==", 4, "Read four residues");
    if(m->num_residues != 4)
        BAIL_OUT("Didn't read any residues: can't complete tests");

    is(m->residues[0].name, "GLU", "Residue 0 is GLU");
    is(m->residues[1].name, "VAL", "Residue 1 is VAL");
    is(m->residues[2].name, "TYR", "Residue 2 is TYR");
    is(m->residues[3].name, "LEU", "Residue 3 is LEU");

    ok(m->linear_springs[0].a->id == m->residues[0].atoms[0].id, "Spring 1 attached to residue 1");
    ok(m->linear_springs[0].b->id == m->residues[1].atoms[0].id, "Spring 1 attached to residue 2");
    ok(abs(m->linear_springs[0].distance - 4.0) < 1e9, "Spring 1 distance is 4.0A");
    ok(abs(m->linear_springs[0].constant - DEFAULT_SPRING_CONSTANT) < 1e9, "Spring 1 constant is %fA", DEFAULT_SPRING_CONSTANT);
    ok(m->linear_springs[0].cutoff < 0, "Spring 1 cutoff disabled");

    ok(m->linear_springs[1].a->id == m->residues[0].atoms[0].id, "Spring 2 attached to residue 1");
    ok(m->linear_springs[1].b->id == m->residues[2].atoms[0].id, "Spring 2 attached to residue 3");
    ok(abs(m->linear_springs[1].distance - 3.0) < 1e9, "Spring 2 distance is 3.0A");
    ok(abs(m->linear_springs[1].constant - 0.2) < 1e9, "Spring 2 constant is 0.2A");
    ok(m->linear_springs[1].cutoff < 0, "Spring 2 cutoff disabled");

    ok(m->linear_springs[2].a->id == m->residues[2].atoms[0].id, "Spring 3 attached to residue 3");
    ok(m->linear_springs[2].b->id == m->residues[3].atoms[0].id, "Spring 3 attached to residue 3");
    ok(abs(m->linear_springs[2].distance - 4.0) < 1e9, "Spring 3 distance is 4.0A");
    ok(abs(m->linear_springs[2].constant - 0.3) < 1e9, "Spring 3 constant is 0.3A");
    ok(abs(m->linear_springs[2].cutoff - 0.1) < 1e9, "Spring 3 cutoff is 0.1A");

    ok(m->torsion_springs[0].a1->id == m->residues[0].atoms[0].id, "Torsion a1 = a1");
    ok(m->torsion_springs[0].a2->id == m->residues[1].atoms[0].id, "Torsion a2 = a2");
    ok(m->torsion_springs[0].a3->id == m->residues[2].atoms[0].id, "Torsion a3 = a3");
    ok(m->torsion_springs[0].a4->id == m->residues[3].atoms[1].id, "Torsion a4 = a4");
    ok(abs(m->torsion_springs[0].angle - 40) < 1e9, "Angle correct");
    ok(abs(m->torsion_springs[0].constant - 0.1) < 1e9, "Constant correct");
}

int main(int argc, char **argv){
    plan(62);

    struct model *ms = springreader_parse_str(pdb);
    check_model(ms);
    model_free(ms);

    //Write temporary file with the same string
    char tmpfile_name[13];
    strcpy(tmpfile_name, "springXXXXXX");
    int tmpfile_fd = mkstemp(tmpfile_name);
    write(tmpfile_fd, pdb, strlen(pdb));

    struct model *mf = springreader_parse_file(tmpfile_name);
    check_model(mf);
    model_free(mf);
    unlink(tmpfile_name);

    done_testing();
}
