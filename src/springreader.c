#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "springreader.h"
#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"

enum section {
    PREAMBLE,
    LINEAR_SPRINGS,
    TORSION_SPRINGS,
    POSITIONS,
    UNKNOWN
};

static enum section parse_section_header(
        const char *line);
static void parse_position_line(
        const char *line, struct model *m);
static void parse_preamble_line(
        const char *line, struct model *m);
static void parse_linear_spring_line(
        const char *line, struct model *m);
static void parse_torsion_spring_line(
        const char *line, struct model *m);
void parse_sequence(const char *seq, struct model *m);
static void parse_line(
        const char *line, struct model *m, enum section *section);
static void scan_double(
        const char *line, const char *value, double *dst);
static bool scan_bool(const char *line, const char *value);

/**
 * Parse  a configuration string into a model.
 *
 * \return A populated model or NULL on error.
 */
struct model * springreader_parse_str(const char *str){
    struct model *m = NULL;
    char *copy = strdup(str);
    if(!copy) goto exit;

    m = model_alloc();
    if(!m) goto free_copy;

    enum section section = PREAMBLE;
    char *line;
    for(line = strtok(copy, "\n"); line; line = strtok(NULL, "\n")){
        parse_line(line, m, &section);
    }
free_copy:
    free(copy);
exit:
    return m;
}

struct model * springreader_parse_file(const char *file){
    struct model *m = model_alloc();
    if(!m) goto bail;

    FILE *f = fopen(file, "r");
    if(!f){
        perror("Error reading input file");
        goto free_model;
    }

    enum section section = PREAMBLE;
    size_t line_sz;
    char *line_buffer = NULL;
    while(getline(&line_buffer, &line_sz, f) != -1){
        parse_line(line_buffer, m, &section);
        free(line_buffer);
        line_buffer = NULL;
    }
    free(line_buffer);

    fclose(f);
    return m;
free_model:
    free(m);
bail:
    return NULL;
}

enum section parse_section_header(const char *line){
    int len = (strrchr(line, ']') - line) - 1;

    char section_name[len+1];
    strncpy(section_name, line+1, len);
    section_name[len] = '\0';

    for(char *p = section_name; *p; p++)
        *p = tolower(*p);

    if(strcmp(section_name, "linear") == 0)
        return LINEAR_SPRINGS;
    else if(strcmp(section_name, "torsion") == 0)
        return TORSION_SPRINGS;
    else if(strcmp(section_name, "position") == 0)
        return POSITIONS;
    return UNKNOWN;
}

void parse_position_line(const char *line, struct model *m){
    int i;
    double x, y, z;
    int n = sscanf(line, "%d %lf %lf %lf", &i, &x, &y, &z);

    if(n != 4){
        fprintf(stderr, "I didn't understand the line `%s'", line);
        return;
    }
    if(i < 1 || (unsigned int) i > m->num_residues){
        fprintf(stderr, "%i is not a valid residue\n", i);
        return;
    }

    i--;
    m->residues[i].atoms[0].position.c[0] = x;
    m->residues[i].atoms[0].position.c[1] = y;
    m->residues[i].atoms[0].position.c[2] = z;
    m->residues[i].atoms[0].synthesised = true;
}

void parse_preamble_line(const char *line, struct model *m){
    char param[strlen(line)];
    char value[strlen(line)];
    sscanf(line, "%s = %s", param, value);

    //Do case-insensitive comparisons of key name
    for(char *p = param; *p; p++) *p = tolower(*p);

    if(strcmp(param, "sequence") == 0){
        parse_sequence(value, m);
    }else if(strcmp(param, "timestep") == 0){
        scan_double(line, value, &m->timestep);
    }else if(strcmp(param, "synth_time") == 0){
        scan_double(line, value, &m->synth_time);
    }else if(strcmp(param, "drag_coefficient") == 0){
        scan_double(line, value, &m->drag_coefficient);
    }else if(strcmp(param, "max_synth_angle") == 0){
        scan_double(line, value, &m->max_synth_angle);
    }else if(strcmp(param, "use_sterics") == 0){
        m->use_sterics = scan_bool(line, value);
    }else if(strcmp(param, "fix") == 0){
        m->fix = scan_bool(line, value);
    }
}

void scan_double(const char *line, const char *value, double *dst){
    if(sscanf(value, "%lf", dst) != 1)
        fprintf(stderr, "Couldn't interpret line %s\n", line);
}

bool scan_bool(const char *line, const char *value){
    if(strcmp(value, "true") == 0){
        return true;
    }else if(strcmp(value, "false") == 0){
        return false;
    }else{
        fprintf(
                stderr,
                "Couldn't interpret boolean `%s' (assuming false): %s",
                value, line);
        return false;
    }
}

void parse_torsion_spring_line(const char *line, struct model *m){
    int r1, r2, r3, r4;
    double angle, constant, cutoff;
    int num_matched = sscanf(line, "%d %d %d %d %lf %lf %lf",
            &r1, &r2, &r3, &r4,
            &angle, &constant, &cutoff);
    switch(num_matched){
        case 6:
            cutoff = -1;
        case 7:
            break;
        default:
            fprintf(stderr, "Couldn't interpret torsion spring: %s\n", line);
            return;
    }

    int r[4] = {r1, r2, r3, r4};
    for(size_t i=0; i < 4; i++){
        if(r[i] <= 0 || r[i] > m->num_residues){
            fprintf(stderr, "Residue %d does not exist -- ignoring spring.\n",
                    r[i]);
            return;
        }
    }

    //Alter the indexing; the file uses 1-indexing, we use 0-indexing
    r1--;
    r2--;
    r3--;
    r4--;

    m->num_torsion_springs++;
    m->torsion_springs = realloc(
            m->torsion_springs,
            m->num_torsion_springs
            * sizeof(*m->torsion_springs));

    torsion_spring_init(&m->torsion_springs[m->num_torsion_springs-1],
            &m->residues[r1].atoms[0], &m->residues[r2].atoms[0],
            &m->residues[r3].atoms[0], &m->residues[r4].atoms[0],
            angle, constant);

}

void parse_linear_spring_line(const char *line, struct model *m){
    int i, j;
    double distance, constant, cutoff;
    int num_matched = sscanf(line, "%d %d %lf %lf %lf",
            &i, &j,
            &distance, &constant, &cutoff);
    switch(num_matched){
        case 3:
            constant = DEFAULT_SPRING_CONSTANT;
        case 4:
            cutoff = -1;
        case 5:
            break;
        default:
            fprintf(stderr, "Couldn't interpret linear spring: %s\n", line);
            return;
    }

    if(i < 1 || i > m->num_residues){
        fprintf(stderr, "Residue %d does not exist -- ignoring spring.\n", i);
        return;
    }else if(j < 1 || i > m->num_residues){
        fprintf(stderr, "Residue %d does not exist -- ignoring spring.\n", j);
        return;
    }
    //Springs are indexed by 1 in the file and zero in the program
    i--;
    j--;

    m->num_linear_springs++;
    m->linear_springs = realloc(
            m->linear_springs,
            m->num_linear_springs
            * sizeof(*m->linear_springs));
    linear_spring_init(
            &m->linear_springs[m->num_linear_springs-1],
            distance, constant,
            &m->residues[i].atoms[0], &m->residues[j].atoms[0]);
}

void parse_line(const char *line, struct model *m, enum section *section){
    if(strlen(line) == 0 || line[0] == '#')
        return;
    else if(line[0] == '['){
        *section = parse_section_header(line);
        if(*section == UNKNOWN)
            fprintf(stderr, "I don't understand the section '%s'\n", line);
    }else{
        switch(*section){
            case PREAMBLE:
                parse_preamble_line(line, m);
                break;
            case LINEAR_SPRINGS:
                parse_linear_spring_line(line, m);
                break;
            case TORSION_SPRINGS:
                parse_torsion_spring_line(line, m);
                break;
            case POSITIONS:
                parse_position_line(line, m);
                break;
            case UNKNOWN:
                //Ignore lines in unknown section
                break;
        }
    }
}

void parse_sequence(const char *seq, struct model *m){

    struct residue *res = malloc(strlen(seq) * sizeof(struct residue));

    //Count sidechains (need to do this because GLY doesn't have one)
    size_t k = 0;
    size_t num_sc = 0;
    size_t i;
    for(i = 0; seq[i]; i++){
        struct residue *r = &res[i];
        struct AA *aa = AA_lookup(&seq[i], 1);
        residue_init(r, i+1);
        strncpy(r->name, aa->threeletter, 3);

        size_t num_atoms = (aa->has_sidechain) ? 2 : 1;
        if(aa->has_sidechain)
            num_sc++;

        r->atoms = malloc(num_atoms * sizeof(struct atom));
        r->num_atoms = num_atoms;
        atom_init(&r->atoms[0], ++k, "CA");
        r->atoms[0].radius = CA_STERIC_RADIUS;
        if(aa->has_sidechain){
            atom_init(&r->atoms[1], ++k, "CB");
            atom_set_AA(&r->atoms[1], aa);
        }
    }

    //Need a spring to attach each CA to each sidechain and to attach each CA
    //to each previous CA.
    struct linear_spring *ls = malloc((num_sc + i) * sizeof(struct linear_spring));
    k = 0;
    for(i=0; seq[i]; i++){
        struct residue *r = &res[i];
        struct AA *aa = AA_lookup(&seq[i], 1);

        if(r->num_atoms > 1){

            for(size_t j=1; j < r->num_atoms; j++){
                linear_spring_init(
                        &ls[k++],
                        aa->sc_bond_len, SC_BB_SPRING_CONSTANT,
                        &r->atoms[0], &r->atoms[j]);
            }
        }
        if(i > 0){
            linear_spring_init(&ls[k++],
                    CA_CA_LEN, BB_BB_SPRING_CONSTANT,
                    &res[i-1].atoms[0], &res[i].atoms[0]);
        }
    }
    m->residues = res;
    m->linear_springs = ls;
    m->num_residues = i;
    m->num_linear_springs = k;
}
