#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "springreader.h"
#include "model.h"
#include "rk4.h"

static struct option opts[] = {
    {"help",     no_argument,       0, 'h'},
    {"snapshot", required_argument, 0, 's'},
    {"until",    required_argument, 0, 'u'},
    {0, 0, 0, 0}
};
const char *opt_str = "hs:u:";

const char *usage_str =
"Usage: poing [OPTIONS] <SPEC>\n"
"\n"
"Argument SPEC is mandatatory and must be a specification file for a "
"simulation.\n"
"Available options:\n"
"  -h, --help         Display this help message.\n"
"  -s, --snapshot=N   Write a PDB snapshot every N steps.\n"
"  -u, --until=T      Run until time T.\n"
;
int snapshot = -1;
int until = 100;

void usage(const char *msg, int exitval){
    FILE *out = (exitval < 2) ? stdout : stderr;
    if(msg)
        fprintf(out, "%s\n", msg);
    fprintf(out, usage_str);
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
            case 'u':
                until = atof(optarg);
                break;
        }
    }
    if(optind >= argc)
        usage("No specification file supplied.", 2);
    return argv[optind];
}

int main(int argc, char **argv){
    char * spec = get_options(argc, argv);
    struct model *model = springreader_parse_file(spec);
    if(!model)
        return 2;

    struct model state;
    while(model->time < until){
        model_synth(&state, model);
        rk4_push(&state);
        if(snapshot > 0 &&
                (int)(state.time / snapshot) > (int)(model->time / snapshot)){
            char *pdb = model_pdb(&state, true);
            printf("MODEL     % d\n", (int)state.time / snapshot);
            printf("%s", pdb);
            printf("ENDMDL\n");
        }
        model->time = state.time;
    }
    return 0;
}
