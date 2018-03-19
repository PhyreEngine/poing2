#include "debug.h"
#include "model.h"
#include "linear_spring.h"
#include "bond_angle.h"
#include "torsion_spring.h"

static bool do_print(struct model *m);

#define DEBUG_LINEAR_FIELDS \
    "time " \
    "a.id a.name b.id b.name " \
    "enabled " \
    "equilibrium distance " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z\n"
#define DEBUG_LINEAR_FMT \
    "%f " \
    "%d %s %d %s " \
    "%s " \
    "%f %f " \
    "%f %f %f " \
    "%f %f %f\n"

#define DEBUG_ANGLE_FIELDS \
    "time " \
    "a.id a.name b.id b.name c.id c.name " \
    "enabled " \
    "equilibrium angle " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z " \
    "force.c.x force.c.y force.c.z\n"
#define DEBUG_ANGLE_FMT \
    "%f " \
    "%d %s %d %s %d %s " \
    "%s " \
    "%f %f "\
    "%f %f %f %f %f %f %f %f %f\n"

#define DEBUG_TORSION_FIELDS \
    "time " \
    "a.id a.name b.id b.name c.id c.name d.id d.name " \
    "enabled " \
    "equilibrium angle " \
    "force.a.x force.a.y force.a.z " \
    "force.b.x force.b.y force.b.z " \
    "force.c.x force.c.y force.c.z " \
    "force.d.x force.d.y force.d.z\n"
#define DEBUG_TORSION_FMT \
    "%f " \
    "%d %s %d %s %d %s %d %s " \
    "%s " \
    "%f %f "\
    "%f %f %f %f %f %f %f %f %f %f %f %f\n"

bool do_print(struct model *m){
    //Only debug if we have just passed the interval
    int nsteps      = (int)(m->time / m->debug->interval);
    int nsteps_prev = (int)((m->time - m->timestep) / m->debug->interval);
    if(nsteps_prev == nsteps)
        return false;
    return true;
}

void debug_begin(struct model *m){
    if(m->debug->linear)
        fprintf(m->debug->linear, DEBUG_LINEAR_FIELDS);

    if(m->debug->angle)
        fprintf(m->debug->angle, DEBUG_ANGLE_FIELDS);

    if(m->debug->torsion)
        fprintf(m->debug->torsion, DEBUG_TORSION_FIELDS);
}

void debug_linear(struct model *m, struct linear_spring *s){
    if(!m->debug->linear || !do_print(m))
        return;

    struct vector force_a, force_b;
    struct vector displacement;
    linear_spring_force(&force_a, &force_b, s);
    vsub(&displacement, &s->a->position, &s->b->position);

    fprintf(m->debug->linear, DEBUG_LINEAR_FMT,
            m->time,
            s->a->id, s->a->name,
            s->b->id, s->b->name,
            (s->enabled ? "enabled" : "disabled"),
            s->distance, vmag(&displacement), 
            force_a.c[0], force_a.c[1], force_a.c[2],
            force_b.c[0], force_b.c[1], force_b.c[2]);
}

void debug_torsion(struct model *m, struct torsion_spring *s){
    if(!m->debug->torsion || !do_print(m))
        return;

    struct vector spring_forces[4];
    torsion_spring_force_new(
            &spring_forces[0],
            &spring_forces[1],
            &spring_forces[2],
            &spring_forces[3],
            s);
    fprintf(m->debug->torsion, DEBUG_TORSION_FMT,
            m->time,
            s->a1->id, s->a1->name,
            s->a2->id, s->a2->name,
            s->a3->id, s->a3->name,
            s->a4->id, s->a4->name,
            (s->enabled ? "enabled" : "disabled"),
            s->angle, torsion_spring_angle(s),
            spring_forces[0].c[0],
            spring_forces[0].c[1],
            spring_forces[0].c[2],
            spring_forces[1].c[0],
            spring_forces[1].c[1],
            spring_forces[1].c[2],
            spring_forces[2].c[0],
            spring_forces[2].c[1],
            spring_forces[2].c[2],
            spring_forces[3].c[0],
            spring_forces[3].c[1],
            spring_forces[3].c[2]);
}

void debug_angle(struct model *m, struct bond_angle_spring *s){
    if(!m->debug->angle || !do_print(m))
        return;

    struct vector spring_forces[3];
    bond_angle_force(
            &spring_forces[0],
            &spring_forces[1],
            &spring_forces[2],
            s);
    fprintf(m->debug->angle, DEBUG_ANGLE_FMT,
            m->time,
            s->a1->id, s->a1->name,
            s->a2->id, s->a2->name,
            s->a3->id, s->a3->name,
            (s->enabled ? "enabled" : "disabled"),
            s->angle, bond_angle_angle(s),
            spring_forces[0].c[0],
            spring_forces[0].c[1],
            spring_forces[0].c[2],
            spring_forces[1].c[0],
            spring_forces[1].c[1],
            spring_forces[1].c[2],
            spring_forces[2].c[0],
            spring_forces[2].c[1],
            spring_forces[2].c[2]);
}


