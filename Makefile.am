#bin_PROGRAMS=a.out
#a_out_SOURCES=src/pppoing.cpp
check_PROGRAMS=test_springreader test_vector \
			   test_linear_spring test_torsion_spring
TESTS=test_springreader test_vector \
	  test_linear_spring test_torsion_spring
BUILT_SOURCES=data/AA.c
SUFFIXES=.gperf

test_springreader_SOURCES=t/springreader.c t/tap.c src/springreader.c \
						src/model.c src/residue.c \
						data/AA.c src/linear_spring.c src/torsion_spring.c \
						src/vector.c

test_vector_SOURCES=t/vector.c src/vector.c t/tap.c

test_linear_spring_SOURCES=t/linear_spring.c src/linear_spring.c \
						   src/vector.c src/residue.c data/AA.c t/tap.c

test_torsion_spring_SOURCES=t/torsion_spring.c src/torsion_spring.c \
						   src/vector.c src/residue.c data/AA.c t/tap.c

.gperf.c:
	$(GPERF) $< --output-file $@
