#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "springreader.h"
#include "residue.h"
#include "linear_spring.h"
#include "torsion_spring.h"
#include "bond_angle.h"
#include "rama.h"

enum section {
    PREAMBLE,
    LINEAR_SPRINGS,
    TORSION_SPRINGS,
    PDB,
    ANGLES,
    RAMA_DATA,
    RAMA,
    UNKNOWN
};

static enum section parse_section_header(
        const char *line);
static void parse_preamble_line(
        const char *line, struct model *m);
static void parse_linear_spring_line(
        const char *line, struct model *m);
static void parse_torsion_spring_line(
        const char *line, struct model *m);
static void parse_bond_angle_line(
        const char *line, struct model *m);
static void parse_rama_data_line(
        const char *line, struct model *m);
static void parse_rama_line(
        const char *line, struct model *m);
static void parse_pdb_line(
        const char *line, struct model *m);
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
    else if(strcmp(section_name, "pdb") == 0)
        return PDB;
    else if(strcmp(section_name, "angle") == 0)
        return ANGLES;
    else if(strcmp(section_name, "ramachandran data") == 0)
        return RAMA_DATA;
    else if(strcmp(section_name, "ramachandran") == 0)
        return RAMA;
    return UNKNOWN;
}

void parse_pdb_line(const char *line, struct model *m){
    size_t len = strlen(line);
    char buf[10];
    buf[0] = '\0';

    if(len < 54)
        return;
    if(strncmp(line, "ATOM", 4) != 0)
        return;

    char atom_name[5];
    char res_name[4];
    size_t atom_id, res_id;
    char chain_id;
    double x, y, z;

    strncat(buf, line + 6, 5);
    sscanf(buf, "%lu", &atom_id);
    buf[0] = '\0';

    strncat(buf, line + 12, 4);
    sscanf(buf, " %s ", atom_name);
    buf[0] = '\0';

    strncat(buf, line + 17, 3);
    sscanf(buf, " %s ", res_name);
    buf[0] = '\0';

    strncat(buf, line + 21, 1);
    sscanf(buf, "%c", &chain_id);
    buf[0] = '\0';

    strncat(buf, line + 22, 4);
    sscanf(buf, "%lu", &res_id);
    buf[0] = '\0';

    strncat(buf, line + 30, 8);
    sscanf(buf, "%lf", &x);
    buf[0] = '\0';

    strncat(buf, line + 38, 8);
    sscanf(buf, "%lf", &y);
    buf[0] = '\0';

    strncat(buf, line + 46, 8);
    sscanf(buf, "%lf", &z);
    buf[0] = '\0';

    if(res_id > m->num_residues){
        struct residue *r = model_push_residue(m, res_id);
        if(!m->do_synthesis)
            r->synthesised = true;
        strcpy(r->name, res_name);
    }
    struct residue *res = &m->residues[m->num_residues-1];

    struct atom *atom = residue_push_atom(res, atom_id, atom_name);
    if(!m->do_synthesis)
        atom->synthesised = true;
    atom->position.c[0] = x;
    atom->position.c[1] = y;
    atom->position.c[2] = z;

    const struct atom_description *desc = atom_description_lookup(
            atom_name, strlen(atom_name));
    atom_set_atom_description(atom, desc);
}

void parse_preamble_line(const char *line, struct model *m){
    char param[strlen(line)];
    char value[strlen(line)];
    sscanf(line, "%s = %s", param, value);

    //Do case-insensitive comparisons of key name
    for(char *p = param; *p; p++) *p = tolower(*p);

    if(strcmp(param, "timestep") == 0){
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
    }else if(strcmp(param, "threestate") == 0){
        m->threestate = scan_bool(line, value);
    }else if(strcmp(param, "use_water") == 0){
        m->use_water = scan_bool(line, value);
    }else if(strcmp(param, "shield_drag") == 0){
        m->shield_drag = scan_bool(line, value);
    }else if(strcmp(param, "do_synthesis") == 0){
        m->do_synthesis = scan_bool(line, value);
    }else if(strcmp(param, "until") == 0){
        scan_double(line, value, &m->until);
    }else{
        fprintf(stderr, "I don't understand the parameter: %s\n", param);
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

void parse_bond_angle_line(const char *line, struct model *m){
    //Residue numbers
    int r1, r2, r3;
    //Atom names
    char na1[4], na2[4], na3[4];
    double angle, constant, cutoff;

    int num_matched = sscanf(line, "%d %s %d %s %d %s %lf %lf %lf",
            &r1, na1, &r2, na2, &r3, na3,
            &angle, &constant, &cutoff);
    switch(num_matched){
        case 8:
            cutoff = -1;
        case 9:
            break;
        default:
            fprintf(stderr, "Couldn't interpret torsion spring: %s\n", line);
            return;
    }

    int r[3] = {r1, r2, r3};
    for(size_t i=0; i < 3; i++){
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

    m->num_bond_angles++;
    m->bond_angles = realloc(
            m->bond_angles,
            m->num_bond_angles * sizeof(*m->bond_angles));
    struct atom *a1, *a2, *a3;
    a1 = a2 = a3 = NULL;

    for(size_t i=0; i < m->residues[r1].num_atoms; i++)
        if(strcmp(m->residues[r1].atoms[i].name, na1) == 0)
            a1 = &m->residues[r1].atoms[i];

    for(size_t i=0; i < m->residues[r2].num_atoms; i++)
        if(strcmp(m->residues[r2].atoms[i].name, na2) == 0)
            a2 = &m->residues[r2].atoms[i];

    for(size_t i=0; i < m->residues[r3].num_atoms; i++)
        if(strcmp(m->residues[r3].atoms[i].name, na3) == 0)
            a3 = &m->residues[r3].atoms[i];

    if(a1 == NULL || a2 == NULL || a3 == NULL){
        fprintf(stderr,
                "Error reading line. "
                "I'm continuing, but this is probably very bad. "
                "Fix your springs file. Offending line:\n"
                "%s",
                line);
    }else{
        bond_angle_spring_init(&m->bond_angles[m->num_bond_angles-1],
                a1, a2, a3,
                angle, constant);
    }

}

static void parse_rama_data_line(
        const char *line, struct model *m){
    //Type
    char type_str[strlen(line)];

    //Check for the rama data
    char filename[strlen(line)];
    if(sscanf(line, "%s = %s", type_str, filename) == 2){
        rama_read_closest(filename, rama_parse_type(type_str));
    }else{
        fprintf(stderr, "Couldn't interpret line: %s\n", line);
    }
}

static void parse_rama_line(
        const char *line, struct model *m){
    //Residue number
    int res_num;
    //Type
    char type_str[strlen(line)];
    float constant = DEFAULT_RAMA_CONST;

    int num_matched = sscanf(line, "%d %s %f", &res_num, type_str, &constant);
    if(num_matched < 2 || num_matched > 3){
        fprintf(stderr, "Couldn't interpret Rama constraint: %s\n", line);
        return;
    }

    if(!rama_is_inited(rama_parse_type(type_str))){
        fprintf(stderr,
                "Ramachandran type %s is not initialised! "
                "Use the 'ramachandran data' section",
                type_str);
        exit(1);
    }

    //For residue i, the phi angle needs the previous residue to exist (res_num
    //>= 2) and the psi angle needs the next residue to exist (res_num <
    //m->num_residues).
    if(res_num < m->num_residues && res_num >= 2){

        m->num_rama_constraints++;
        m->rama_constraints = realloc(
                m->rama_constraints,
                m->num_rama_constraints
                * sizeof(*m->rama_constraints));

        //res_num is 1-indexed, but we used 0-indexing
        struct residue *res      = &m->residues[res_num - 1];
        struct residue *next_res = &m->residues[res_num];
        struct residue *prev_res = &m->residues[res_num - 2];

        rama_init(
                &m->rama_constraints[m->num_rama_constraints - 1],
                res, next_res, prev_res, type_str, constant);
        rama_random_init(&m->rama_constraints[m->num_rama_constraints - 1]);
    }
}

void parse_torsion_spring_line(const char *line, struct model *m){
    //Residue numbers
    int r1, r2, r3, r4;
    //Atom names
    char na1[4], na2[4], na3[4], na4[4];
    double angle, constant, cutoff;

    int num_matched = sscanf(line, "%d %s %d %s %d %s %d %s %lf %lf %lf",
            &r1, na1, &r2, na2, &r3, na3, &r4, na4,
            &angle, &constant, &cutoff);
    switch(num_matched){
        case 10:
            cutoff = -1;
        case 11:
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
    struct atom *a1, *a2, *a3, *a4;
    a1 = a2 = a3 = a4 = NULL;

    for(size_t i=0; i < m->residues[r1].num_atoms; i++)
        if(strcmp(m->residues[r1].atoms[i].name, na1) == 0)
            a1 = &m->residues[r1].atoms[i];

    for(size_t i=0; i < m->residues[r2].num_atoms; i++)
        if(strcmp(m->residues[r2].atoms[i].name, na2) == 0)
            a2 = &m->residues[r2].atoms[i];

    for(size_t i=0; i < m->residues[r3].num_atoms; i++)
        if(strcmp(m->residues[r3].atoms[i].name, na3) == 0)
            a3 = &m->residues[r3].atoms[i];

    for(size_t i=0; i < m->residues[r4].num_atoms; i++)
        if(strcmp(m->residues[r4].atoms[i].name, na4) == 0)
            a4 = &m->residues[r4].atoms[i];

    if(a1 == NULL || a2 == NULL || a3 == NULL || a4 == NULL){
        fprintf(stderr,
                "Error reading line. "
                "I'm continuing, but this is probably very bad. "
                "Fix your springs file. Offending line:\n"
                "%s",
                line);
    }else{
        torsion_spring_init(&m->torsion_springs[m->num_torsion_springs-1],
                a1, a2, a3, a4,
                angle, constant);
    }

}

void parse_linear_spring_line(const char *line, struct model *m){
    int i, j;
    char i_name[5], j_name[5];
    double distance, constant, cutoff;
    int num_matched = sscanf(line, "%d %4s %d %4s %lf %lf %lf",
            &i, i_name,
            &j, j_name,
            &distance, &constant, &cutoff);
    switch(num_matched){
        case 5:
            constant = DEFAULT_SPRING_CONSTANT;
        case 6:
            cutoff = -1;
        case 7:
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

    struct atom *atom_i, *atom_j;
    for(int k=0; k < m->residues[i].num_atoms; k++){
        if(strcmp(i_name, m->residues[i].atoms[k].name) == 0){
            atom_i = &m->residues[i].atoms[k];
        }
    }
    for(int k=0; k < m->residues[j].num_atoms; k++){
        if(strcmp(j_name, m->residues[j].atoms[k].name) == 0){
            atom_j = &m->residues[j].atoms[k];
        }
    }

    m->num_linear_springs++;
    m->linear_springs = realloc(
            m->linear_springs,
            m->num_linear_springs
            * sizeof(*m->linear_springs));
    linear_spring_init(
            &m->linear_springs[m->num_linear_springs-1],
            distance, constant,
            atom_i, atom_j);
    m->linear_springs[m->num_linear_springs-1].cutoff = cutoff;
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
            case PDB:
                parse_pdb_line(line, m);
                break;
            case ANGLES:
                parse_bond_angle_line(line, m);
                break;
            case RAMA_DATA:
                parse_rama_data_line(line, m);
                break;
            case RAMA:
                parse_rama_line(line, m);
                break;
            case UNKNOWN:
                //Ignore lines in unknown section
                break;
        }
    }
}
