#ifndef PROFILE_H_
#define PROFILE_H_

#include <time.h>
#include <stdio.h>

struct profile {
    FILE *out;
    struct timespec start;
};

void profile_init(struct profile *profiler, FILE *out);
void profile_start(struct profile *profiler);
long long profile_duration(struct profile *profiler);
void profile_end_internal(struct profile *profiler, const char *fmt, ...);

//Macro required for simple re-ordering of va_list parameters. Required because
//we want to pass the time parameter to the format string as the last item.
#define profile_end(profiler, fmt, ...)     \
    profile_end_internal(                   \
        (profiler),                         \
        (fmt), __VA_ARGS__,                 \
        profile_duration((profiler)))

#endif /* PROFILE_H_ */
