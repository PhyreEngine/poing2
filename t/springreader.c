#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/springreader.h"
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/linear_spring.h"
#include "../src/torsion_spring.h"
#include "../src/bond_angle.h"
#include "tap.h"

const char *json =
"{\n"
"    \"sequence\": \"EVYL\",\n"
"    \"timestep\": 0.01,\n"
"    \"synth_time\": 100,\n"
"    \"drag_coefficient\": 0.1,\n"
"    \"atoms\": [\n"
"        {\"id\": 1, \"name\": \"CA\",  \"residue\": 1},\n"
"        {\"id\": 2, \"name\": \"GLU\", \"residue\": 1},\n"
"        {\"id\": 3, \"name\": \"CA\",  \"residue\": 2},\n"
"        {\"id\": 4, \"name\": \"VAL\", \"residue\": 2},\n"
"        {\"id\": 5, \"name\": \"CA\",  \"residue\": 3},\n"
"        {\"id\": 6, \"name\": \"TYR\", \"residue\": 3},\n"
"        {\"id\": 7, \"name\": \"CA\",  \"residue\": 4},\n"
"        {\"id\": 8, \"name\": \"LEU\", \"residue\": 4}\n"
"    ],\n"
"    \"linear\": [\n"
"        {\n"
"            \"atoms\": [1, 3],\n"
"            \"distance\": 4.0\n"
"        },\n"
"        {\n"
"            \"atoms\": [1, 5],\n"
"            \"distance\": 3.0,\n"
"            \"constant\": 0.2\n"
"        },\n"
"        {\n"
"            \"atoms\": [5, 7],\n"
"            \"distance\": 4.0,\n"
"            \"constant\": 0.3,\n"
"            \"cutoff\": 0.1\n"
"        }\n"
"    ],\n"
"    \"torsion\": [\n"
"        {\n"
"            \"atoms\": [1, 3, 5, 8],\n"
"            \"angle\": 40,\n"
"            \"constant\": 0.1\n"
"        }\n"
"    ],\n"
"    \"angle\": [\n"
"        {\n"
"            \"atoms\": [1, 3, 5],\n"
"            \"angle\": 45,\n"
"            \"constant\": 0.1\n"
"        }\n"
"    ],\n"
"    \"ramachandran\": {\n"
"        \"data\": {\n"
"            \"general\": \"data/boundary-general-nosec.data\",\n"
"            \"alpha\": \"data/boundary-alpha.data\",\n"
"            \"beta\": \"data/boundary-beta.data\",\n"
"            \"proline\": \"data/boundary-pro.data\",\n"
"            \"pre_proline\": \"data/boundary-prepro.data\",\n"
"            \"glycine\": \"data/boundary-gly-sym-nosec.data\",\n"
"            \"alanine\": \"data/boundary-ala-nosec.data\"\n"
"        },\n"
"        \"constraints\": [\n"
"            {\"residue\": 2, \"type\": \"general\"},\n"
"            {\"residue\": 3, \"type\": \"alanine\"}\n"
"        ]\n"
"    }\n"
"}\n";

void check_model(struct model *m){
    fis(m->timestep, 0.01, 1e-10, "Timestep");
    fis(m->synth_time, 100, 1e-10, "Synth time");
    fis(m->drag_coefficient, 0.1, 1e-10, "Drag coefficient");

    //3 springs
    cmp_ok(m->num_linear_springs, "==",  3, "Read three linear springs");
    cmp_ok(m->num_torsion_springs, "==", 1, "Read one torsion spring");
    cmp_ok(m->num_residues, "==", 4, "Read four residues");
    cmp_ok(m->num_rama_constraints, "==", 2, "Read two rama constraints");
    if(m->num_residues != 4)
        BAIL_OUT("Didn't read any residues: can't complete tests");

    is(m->residues[0].name, "GLU", "Residue 0 is GLU");
    is(m->residues[1].name, "VAL", "Residue 1 is VAL");
    is(m->residues[2].name, "TYR", "Residue 2 is TYR");
    is(m->residues[3].name, "LEU", "Residue 3 is LEU");

    ok(m->linear_springs[0].a->id == m->atoms[0].id, "Spring 1 attached to residue 1");
    ok(m->linear_springs[0].b->id == m->atoms[2].id, "Spring 1 attached to residue 2");
    ok(abs(m->linear_springs[0].distance - 4.0) < 1e9, "Spring 1 distance is 4.0A");
    ok(abs(m->linear_springs[0].constant - DEFAULT_SPRING_CONSTANT) < 1e9, "Spring 1 constant is %fA", DEFAULT_SPRING_CONSTANT);
    ok(m->linear_springs[0].cutoff < 0, "Spring 1 cutoff disabled");

    ok(m->linear_springs[1].a->id == m->atoms[0].id, "Spring 2 attached to residue 1");
    ok(m->linear_springs[1].b->id == m->atoms[4].id, "Spring 2 attached to residue 3");
    ok(abs(m->linear_springs[1].distance - 3.0) < 1e9, "Spring 2 distance is 3.0A");
    ok(abs(m->linear_springs[1].constant - 0.2) < 1e9, "Spring 2 constant is 0.2A");
    ok(m->linear_springs[1].cutoff < 0, "Spring 2 cutoff disabled");

    ok(m->linear_springs[2].a->id == m->atoms[4].id, "Spring 3 attached to residue 3");
    ok(m->linear_springs[2].b->id == m->atoms[6].id, "Spring 3 attached to residue 3");
    ok(abs(m->linear_springs[2].distance - 4.0) < 1e9, "Spring 3 distance is 4.0A");
    ok(abs(m->linear_springs[2].constant - 0.3) < 1e9, "Spring 3 constant is 0.3A");
    ok(abs(m->linear_springs[2].cutoff - 0.1) < 1e9, "Spring 3 cutoff is 0.1A");

    ok(m->torsion_springs[0].a1->id == m->atoms[0].id, "Torsion a1 = a1");
    ok(m->torsion_springs[0].a2->id == m->atoms[2].id, "Torsion a2 = a2");
    ok(m->torsion_springs[0].a3->id == m->atoms[4].id, "Torsion a3 = a3");
    ok(m->torsion_springs[0].a4->id == m->atoms[7].id, "Torsion a4 = a4");
    ok(abs(m->torsion_springs[0].angle - 40) < 1e9, "Angle correct");
    ok(abs(m->torsion_springs[0].constant - 0.1) < 1e9, "Constant correct");

    ok(m->bond_angles[0].a1->id == m->atoms[0].id, "Angle a1 = a1");
    ok(m->bond_angles[0].a2->id == m->atoms[2].id, "Angle a2 = a2");
    ok(m->bond_angles[0].a3->id == m->atoms[4].id, "Angle a3 = a3");
    ok(abs(m->bond_angles[0].angle - 45) < 1e9, "Angle correct");
    ok(abs(m->bond_angles[0].constant - 0.1) < 1e9, "Constant correct");
}

int main(int argc, char **argv){
    plan(74);

    struct model *ms = springreader_parse_str(json);
    if(!ms)
        BAIL_OUT("Error reading valid JSON string\n");

    check_model(ms);
    model_free(ms);

    //Write temporary file with the same string
    char tmpfile_name[13];
    strcpy(tmpfile_name, "springXXXXXX");
    int tmpfile_fd = mkstemp(tmpfile_name);
    write(tmpfile_fd, json, strlen(json));

    struct model *mf = springreader_parse_file(tmpfile_name);
    if(!mf)
        BAIL_OUT("Error reading valid JSON file\n");

    check_model(mf);
    model_free(mf);
    unlink(tmpfile_name);

    done_testing();
}
