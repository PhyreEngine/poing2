#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/springreader.h"
#include "../src/model.h"
#include "tap.h"

const char *pdb =
"#Sequence:\n"
"sequence = PAPLEQMR\n"
"timestep = 0.01\n"
"synth_time = 100\n"
"drag_coefficient = 0.1\n"
"#List of springs (i, j, distance)\n"
"[Linear]\n"
"1 2 4.0 0.3\n"
"1 3 3.0 0.2\n"
"3 4 4.0 0.3\n"
"[Torsion]\n"
"1 2 3 4 40 0.1\n";

void check_model(struct model *m){
    fis(m->timestep, 0.01, 1e-10, "Timestep");
    fis(m->synth_time, 100, 1e-10, "Synth time");
    fis(m->drag_coefficient, 0.1, 1e-10, "Drag coefficient");

    cmp_ok(m->num_linear_springs, "==", 3, "Read three linear springs");
    cmp_ok(m->num_torsion_springs, "==", 1, "Read one torsion spring");
    cmp_ok(m->num_residues, "==", 8, "Read eight residues");
    if(m->num_residues != 8)
        BAIL_OUT("Didn't read any residues: can't complete tests");

    is(m->residues[0].aa->oneletter, "P", "Residue 0 is P");
    is(m->residues[1].aa->oneletter, "A", "Residue 0 is A");
    is(m->residues[2].aa->oneletter, "P", "Residue 0 is P");
    is(m->residues[3].aa->oneletter, "L", "Residue 0 is L");
    is(m->residues[4].aa->oneletter, "E", "Residue 0 is E");
    is(m->residues[5].aa->oneletter, "Q", "Residue 0 is Q");
    is(m->residues[6].aa->oneletter, "M", "Residue 0 is M");
    is(m->residues[7].aa->oneletter, "R", "Residue 0 is R");

    ok(m->linear_springs[0].a == &m->residues[0], "Spring 1 attached to residue 1");
    ok(m->linear_springs[0].b == &m->residues[1], "Spring 1 attached to residue 2");
    ok(abs(m->linear_springs[0].distance - 4.0) < 1e9, "Spring 1 distance is 4.0A");
    ok(abs(m->linear_springs[0].constant - 0.3) < 1e9, "Spring 1 constant is 0.3A");

    ok(m->linear_springs[1].a == &m->residues[0], "Spring 2 attached to residue 1");
    ok(m->linear_springs[1].b == &m->residues[2], "Spring 2 attached to residue 3");
    ok(abs(m->linear_springs[1].distance - 3.0) < 1e9, "Spring 2 distance is 3.0A");
    ok(abs(m->linear_springs[1].constant - 0.2) < 1e9, "Spring 2 constant is 0.2A");

    ok(m->linear_springs[2].a == &m->residues[2], "Spring 3 attached to residue 3");
    ok(m->linear_springs[2].b == &m->residues[3], "Spring 3 attached to residue 3");
    ok(abs(m->linear_springs[2].distance - 4.0) < 1e9, "Spring 3 distance is 4.0A");
    ok(abs(m->linear_springs[2].constant - 0.3) < 1e9, "Spring 3 constant is 0.3A");

    ok(m->torsion_springs[0].r1 == &m->residues[0], "Torsion r1 = r1");
    ok(m->torsion_springs[0].r2 == &m->residues[1], "Torsion r2 = r2");
    ok(m->torsion_springs[0].r3 == &m->residues[2], "Torsion r3 = r3");
    ok(m->torsion_springs[0].r4 == &m->residues[3], "Torsion r4 = r4");
    ok(abs(m->torsion_springs[0].angle - 40) < 1e9, "Angle correct");
    ok(abs(m->torsion_springs[0].constant - 0.1) < 1e9, "Constant correct");
}

int main(int argc, char **argv){
    plan(64);

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
