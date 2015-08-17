#include <config.h>
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/rk4.h"
#include "../src/springreader.h"
#include "../src/model.h"
#include "../src/residue.h"
#include "../src/linear_spring.h"
#include "../src/torsion_spring.h"
#include "../src/vector.h"

void run(const char *spec, const char *out){
    struct model state;
    struct model *m = springreader_parse_str(spec);

    FILE *fout = fopen(out, "w");
    for(int i=0; i < 10000; i++){
        model_synth(&state, m);
        if(i % 10 == 0){
            fprintf(fout, "MODEL     % d\n", i);
            model_pdb(fout, &state, true);
            fprintf(fout, "ENDMDL\n");
        }
        rk4_push(&state);
        m->time = state.time;
    }
    fclose(fout);
    model_free(m);
}

int main(int argc, char **argv){
#ifdef _GNU_SOURCE
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    const char *lin_spec =
        "timestep = 0.01\n"
        "synth_time = 1\n"
        "shield_drag = false\n"
        "drag_coefficient = -0.5\n"
        "[PDB]\n"
        "ATOM      1  CA  GLY A   1      -1.000   0.000   0.000\n"
        "ATOM      2  CA  GLY A   2       0.000   0.000   0.000\n"
        "ATOM      3  CA  GLY A   3       0.000   0.000   1.000\n"
        "ATOM      4  CA  GLY A   4       0.000   1.000   1.000\n"
        "[Linear]\n"
        "1 CA 2 CA 3.8 1.0\n"
        "2 CA 3 CA 3.8 1.0\n"
        "3 CA 4 CA 3.8 1.0\n"
        ;

    const char *angle_spec =
        "timestep = 0.01\n"
        "synth_time = 1\n"
        "shield_drag = false\n"
        "drag_coefficient = -0.5\n"
        "[PDB]\n"
        "ATOM      1  CA  GLY A   1      -1.000   0.000   0.000\n"
        "ATOM      2  CA  GLY A   2       0.000   0.000   0.000\n"
        "ATOM      3  CA  GLY A   3       0.000   0.000   1.000\n"
        "ATOM      4  CA  GLY A   4       0.000   1.000   1.000\n"
        "[Linear]\n"
        "1 CA 2 CA 3.8 1.0\n"
        "2 CA 3 CA 3.8 1.0\n"
        "3 CA 4 CA 3.8 1.0\n"
        "[Angle]\n"
        "1 CA 2 CA 3 CA 90 0.1\n"
        "2 CA 3 CA 4 CA 90 0.1\n"
        ;

    const char *torsion_spec =
        "timestep = 0.01\n"
        "synth_time = 1\n"
        "shield_drag = false\n"
        "drag_coefficient = -0.5\n"
        "[PDB]\n"
        "ATOM      1  CA  GLY A   1      -1.000   0.000   0.000\n"
        "ATOM      2  CA  GLY A   2       0.000   0.000   0.000\n"
        "ATOM      3  CA  GLY A   3       0.000   0.000   1.000\n"
        "ATOM      4  CA  GLY A   4       0.000   1.000   1.000\n"
        "[Linear]\n"
        "1 CA 2 CA 3.8 1.0\n"
        "2 CA 3 CA 3.8 1.0\n"
        "3 CA 4 CA 3.8 1.0\n"
        "[Angle]\n"
        "1 CA 2 CA 3 CA 90 0.1\n"
        "2 CA 3 CA 4 CA 90 0.1\n"
        "[Torsion]\n"
        "1 CA 2 CA 3 CA 4 CA 45 0.5\n"
        ;

    run(lin_spec, "linear_spring.csv");
    run(angle_spec, "angle_spring.csv");
    run(torsion_spec, "torsion_spring.csv");
}
