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
#include "rattle.h"
#include "sterics.h"
#include "linear_spring.h"

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

    //Make a copy of our model to act as the current state
    struct model state;
    memcpy(&state, model, sizeof(state));

    //Count the number of snapshots written
    size_t num_snapshots = 0;

    //If we are simulating synthesis, the current state will start with no
    //atoms and no residues.
    if(model->do_synthesis){
        state.num_atoms = 0;
        state.num_residues = 0;
    }

    unsigned long i=0;
    while(state.time < state.until){
        //Calculate the number of synthesised atoms from the current time
        int num_synthed = (int)(state.time / state.synth_time) + 1;

        //If we have too few atoms, synthesise the next one
        if(num_synthed > state.num_atoms && state.num_atoms < model->num_atoms){
            size_t new_atom_idx = state.num_atoms;
            state.num_atoms++;
            state.num_residues = state.atoms[new_atom_idx].residue_idx + 1;
            model_synth_atom(&state, new_atom_idx, DEFAULT_MAX_SYNTH_ANGLE);
        }

        //Write PDB file if required
        if(snapshot > 0 && (int)(state.time / snapshot) > num_snapshots){
            printf("MODEL     %lu\n", ++num_snapshots);
            model_pdb(stdout, &state, print_connect);
            printf("ENDMDL\n");
        }

        //Push atoms
        rattle_push(&state);
    }


    steric_grid_free(steric_grid);
    model_free(model);
    return 0;
}
