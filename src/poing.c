#include <config.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "springreader.h"
#include "model.h"
#include "rk4.h"
#include "leapfrog.h"
#include "sterics.h"

#ifdef _GNU_SOURCE
#   include <fenv.h>
#endif

enum state {FROZEN, NORMAL};

static void debug_file(FILE **f, const char *loc, const char *header);

static struct option opts[] = { {"help",     no_argument,       0, 'h'},
    {"snapshot",   required_argument, 0, 's'},
    {"no-connect", no_argument,       0, 'c'},
    {"seed",       required_argument, 0, 'r'},
    {"kinetic",    required_argument, 0, 'k'},
    {"debug-linear",  required_argument, 0, 'l'},
    {"debug-angle",   required_argument, 0, 'a'},
    {"debug-torsion", required_argument, 0, 't'},
    {0, 0, 0, 0}
};
const char *opt_str = "hs:u:s:k:r:l:a:t:";

const char *usage_str =
"Usage: poing [OPTIONS] <SPEC>\n"
"\n"
"Argument SPEC is mandatatory and must be a specification file for a "
"simulation.\n"
"Available options:\n"
"  -h, --help         Display this help message.\n"
"  -s, --snapshot=N   Write a PDB snapshot every N steps.\n"
"  -r, --seed=S       Use fixed random seed S.\n"
"  -k, --kinetic=F    Write kinetic energies to file F.\n"
"      --no-connect   Do not print CONECT records for each spring.\n"
;
int snapshot = -1;
bool print_connect = true;
bool fixed_seed = false;
int random_seed = 0;
char *kinetic = NULL;

struct model_debug debug_opts = {NULL, NULL, NULL, NULL, 0, 0};

void usage(const char *msg, int exitval){
    FILE *out = (exitval < 2) ? stdout : stderr;
    if(msg)
        fprintf(out, "%s\n", msg);
    fprintf(out, "%s", usage_str);
    exit(exitval);
}

char * get_options(int argc, char **argv){
    int c;
    int option_index;
    while((c = getopt_long(argc, argv, opt_str, opts, &option_index)) != -1){
        switch(c){
            case 'h':
                usage(NULL, 1);
            case 's':
                snapshot = atoi(optarg);
                break;
            case 'c':
                print_connect = false;
                break;
            case 'r':
                fixed_seed = true;
                random_seed = atoi(optarg);
                break;
            case 'k':
                kinetic = optarg;
                break;
            case 'l':
                debug_file(&debug_opts.linear, optarg, DEBUG_LINEAR_FIELDS);
                break;
            case 'a':
                debug_file(&debug_opts.angle, optarg, DEBUG_ANGLE_FIELDS);
                break;
            case 't':
                debug_file(&debug_opts.torsion, optarg, DEBUG_TORSION_FIELDS);
                break;
        }
    }
    if(optind >= argc)
        usage("No specification file supplied.", 2);
    return argv[optind];
}

void debug_file(FILE **f, const char *loc, const char *header){
    *f = fopen(loc, "w");
    if(!(*f))
        perror("Error opening debug file");
    else
        fprintf(*f, header);
}

int main(int argc, char **argv){
#if defined(_GNU_SOURCE) && !defined(__FAST_MATH__)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    /* Get options and whatnot */
    char * spec = get_options(argc, argv);
    struct model *model = springreader_parse_file(spec);
    if(!model)
        return 2;
    debug_opts.interval = snapshot;
    model->debug = &debug_opts;

    if(!fixed_seed){
        srand(time(NULL) * getpid());
    }else{
        srand(random_seed);
    }

    /* Open output file for kinetic energy if specified. */
    FILE *kinetic_out = NULL;
    if(kinetic){
        kinetic_out = fopen(kinetic, "w");
        if(!kinetic_out){
            perror("Couldn't open kinetic energy output file");
            exit(1);
        }
    }


    /* Initialise steric grid if it is being used. */
    struct steric_grid *steric_grid = NULL;
    if(model->use_sterics || model->use_water || model->shield_drag){
         steric_grid = steric_grid_alloc(6);
        model->steric_grid = steric_grid;
    }

        /* If we're doing a three-state synthesis then fix/unfix residues
     * and enable/disable springs. */
    enum state three_state = NORMAL;

    /* Run. */
    struct model state;
    unsigned long i=0;
    while(model->time < model->until){

        /* Do synth steps */
        model_synth(&state, model);
        leapfrog_push(&state);


        if(model->threestate && state.num_residues > 0){
            int num_syn  = state.time / state.synth_time;
            double dt    = state.timestep;
            double since = state.time - (num_syn * state.synth_time);
            struct residue *res  = &state.residues[state.num_residues - 1];
            struct residue *prev = NULL;
            if(state.num_residues > 1)
                prev = &state.residues[state.num_residues - 2];

            /*
             * We begin by fixing all atoms except those in the newly
             * synthesised residue, just to make our life a bit easier. Then,
             * we disable all springs except back to the previous residue.
             *
             * After a little bit (1/2 of the synthesis time, we unfreeze the
             * atoms and re-enable springs.
             */
            if(since < model->synth_time / 2 && three_state == NORMAL){
                //Fix atoms
                for(size_t i=0; i < state.num_residues - 1; i++){
                    for(size_t j=0; j < state.residues[i].num_atoms; j++){
                        state.residues[i].atoms[j].fixed = true;
                    }
                }

                //Disable springs
                struct residue *recent = &state.residues[state.num_residues - 1];
                for(size_t i=0; i < state.num_linear_springs; i++){

                    //Allow springs only if they are part of the backbone
                    struct atom *a1 = state.linear_springs[i].a;
                    struct atom *a2 = state.linear_springs[i].b;

                    if(strcmp(a1->name, "CA") == 0 && strcmp(a2->name, "CA") == 0)
                        state.linear_springs[i].enabled = false;
                }
                //Set state
                three_state = FROZEN;
            }else if(since > model->synth_time / 2 && three_state == FROZEN){
                for(size_t i=0; i < state.num_residues - 1; i++){
                    for(size_t j=0; j < state.residues[i].num_atoms; j++){
                        state.residues[i].atoms[j].fixed = false;
                    }
                }

                for(size_t i=0; i < state.num_linear_springs; i++){
                    state.linear_springs[i].enabled = true;
                }
                three_state = NORMAL;
            }
        }



        /* Print snapshot at the correct times. */
        if(snapshot > 0 &&
                (int)(state.time / snapshot) > (int)(model->time / snapshot)){
            printf("MODEL     %lu\n", ++i);
            model_pdb(stdout, &state, print_connect);
            printf("ENDMDL\n");

            /* Also print KE if the output file was supplied. */
            if(kinetic){
                double ke = 0;
                for(size_t j=0; j < state.num_residues; j++){
                    for(size_t k=0; k < state.residues[j].num_atoms; k++){
                        struct atom *atom = &state.residues[j].atoms[k];
                        ke += 0.5* atom->mass * vdot(
                                &atom->velocity, &atom->velocity);
                    }
                }
                fprintf(kinetic_out, "%f\t%f\n", state.time, ke);
                fflush(kinetic_out);
            }
        }
        model->time = state.time;
    }
    if(kinetic)
        fclose(kinetic_out);
    steric_grid_free(steric_grid);
    model_free(model);
    return 0;
}
