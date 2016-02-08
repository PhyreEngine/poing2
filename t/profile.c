#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../src/profile.h"
#include "tap.h"

int main(int argc, char **argv){
    plan(6);

    char buffer[1024];
    FILE *mem = fmemopen(buffer, 1024, "w");
    struct profile profiler;

    profile_init(&profiler, mem);

    profile_start(&profiler);
    sleep(1);
    profile_end(&profiler, "%d\t%s\t%lld\n", 1, "TEST");
    fclose(mem);

    int read_first;
    char read_second[10];
    long long read_third = 0;
    sscanf(buffer, "%d\t%s\t%lld", &read_first, read_second, &read_third);

    cmp_ok(read_first, "==", 1, "Read integer");
    is(read_second, "TEST", "Read string");
    ok(read_third > 0, "Read time");

    done_testing();
}
