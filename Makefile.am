#bin_PROGRAMS=a.out
#a_out_SOURCES=src/pppoing.cpp
check_PROGRAMS=test_springreader
TESTS=test_springreader
BUILT_SOURCES=data/AA.c
SUFFIXES=.gperf

test_springreader_SOURCES=t/springreader.c t/tap.c src/springreader.c \
						src/model.c src/residue.c \
						data/AA.c src/linear_spring.c src/torsion_spring.c

.gperf.c:
	gperf $< --output-file $@
