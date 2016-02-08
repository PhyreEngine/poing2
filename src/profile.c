#include "profile.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

void profile_init(struct profile *profiler, FILE *out){
    profiler->out = out;
}

void profile_start(struct profile *profiler){
    clock_gettime(CLOCK_MONOTONIC_RAW, &profiler->start);
}

long long profile_duration(struct profile *profiler){
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    //Get time
    long long ns_start = profiler->start.tv_sec * 1e9 + profiler->start.tv_nsec;
    long long ns_end   = end.tv_sec * 1e9 + end.tv_nsec;
    return ns_end - ns_start;
}

void profile_end_internal(struct profile *profiler, const char *fmt, ...){

    //Get extra args
    va_list ap;
    va_start(ap, fmt);

    //Print time and extra args
    vfprintf(profiler->out, fmt, ap);

    va_end(ap);

    //Record start again so we don't have to keep calling profile_start
    clock_gettime(CLOCK_MONOTONIC_RAW, &profiler->start);
}
