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
    UNKNOWN
};

static enum section parse_section_header(
        const char *line);
static void parse_preamble_line(
        const char *line, struct model *m);
static struct linear_spring * parse_linear_spring_line(
        const char *line, struct model *m);
static struct torsion_spring * parse_torsion_spring_line(
        const char *line, struct model *m);
static size_t parse_sequence(
        const char *seq, struct model *m, struct residue **res);
static void parse_line(
        const char *line, struct model *m, enum section *section);
static void scan_double(
        const char *line, const char *value, double *dst);

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
    char line_buffer[80];
    while(fgets(line_buffer, 80, f))
        parse_line(line_buffer, m, &section);

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
    return UNKNOWN;
}

void parse_preamble_line(const char *line, struct model *m){
    char param[strlen(line)];
    char value[strlen(line)];
    sscanf(line, "%s = %s", param, value);

    //Do case-insensitive comparisons of key name
    for(char *p = param; *p; p++) *p = tolower(*p);

    if(strcmp(param, "sequence") == 0){
        struct residue *res;
        size_t num_res = parse_sequence(value, m, &res);
        m->num_residues = num_res;
        m->residues = res;
    }else if(strcmp(param, "timestep") == 0){
        scan_double(line, value, &m->timestep);
    }else if(strcmp(param, "synth_time") == 0){
        scan_double(line, value, &m->synth_time);
    }else if(strcmp(param, "drag_coefficient") == 0){
        scan_double(line, value, &m->drag_coefficient);
    }
}

void scan_double(const char *line, const char *value, double *dst){
    if(sscanf(value, "%lf", dst) != 1)
        fprintf(stderr, "Couldn't interpret line %s\n", line);
}

struct torsion_spring * parse_torsion_spring_line(const char *line, struct model *m){
    int r1, r2, r3, r4;
    double angle, constant;
    int num_matched = sscanf(line, "%d %d %d %d %lf %lf",
            &r1, &r2, &r3, &r4,
            &angle, &constant);
    if(num_matched != 6){
        fprintf(stderr, "Couldn't interpret torsion spring: %s\n", line);
        return NULL;
    }
    r1--;
    r2--;
    r3--;
    r4--;
    struct torsion_spring *s = torsion_spring_alloc(
            &m->residues[r1], &m->residues[r2],
            &m->residues[r3], &m->residues[r4],
            angle, constant);
    return s;
}

struct linear_spring * parse_linear_spring_line(const char *line, struct model *m){
    int i, j;
    double distance, constant;
    int num_matched = sscanf(line, "%d %d %lf %lf",
            &i, &j,
            &distance, &constant);
    if(num_matched != 4){
        fprintf(stderr, "Couldn't interpret linear spring: %s\n", line);
        return NULL;
    }
    if(i > m->num_residues || i < 1){
        fprintf(stderr, "Residue %d does not exist -- ignoring spring.\n", i);
        return NULL;
    }else if(j > m->num_residues || j < 1){
        fprintf(stderr, "Residue %d does not exist -- ignoring spring.\n", j);
        return NULL;
    }
    i--;
    j--;

    struct linear_spring *s = linear_spring_alloc(
            distance, constant,
            &m->residues[i], &m->residues[j]
    );
    return s;
}

void parse_line(const char *line, struct model *m, enum section *section){
    struct linear_spring *ls;
    struct torsion_spring *ts;
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
                ls = parse_linear_spring_line(line, m);
                if(ls){
                    m->num_linear_springs++;
                    m->linear_springs = realloc(
                            m->linear_springs,
                            m->num_linear_springs
                            * sizeof(*m->linear_springs));
                    m->linear_springs[m->num_linear_springs-1] = *ls;
                }
                break;
            case TORSION_SPRINGS:
                ts = parse_torsion_spring_line(line, m);
                if(ts){
                    m->num_torsion_springs++;
                    m->torsion_springs = realloc(
                            m->torsion_springs,
                            m->num_torsion_springs
                            * sizeof(*m->torsion_springs));
                    m->torsion_springs[m->num_torsion_springs-1] = *ts;
                }
                break;
            case UNKNOWN:
                //Ignore lines in unknown section
                break;
        }
    }
}

size_t parse_sequence(const char *seq, struct model *m, struct residue **res){
    *res = malloc(strlen(seq) * sizeof(struct residue));
    if(!res)
        return 0;

    size_t i;
    for(i = 0; seq[i]; i++){
        struct AA *aa = AA_lookup(&seq[i], 1);
        residue_init(&(*res)[i], aa, i+1);
    }
    return i;
}
