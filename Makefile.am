bin_PROGRAMS=poing2
poing2_CFLAGS=$(OPENMP_CFLAGS)
poing2_SOURCES=src/poing.c src/rk4.c src/model.c src/residue.c \
			   src/linear_spring.c src/torsion_spring.c src/springreader.c \
			   src/vector.c src/sterics.c data/AA.c


check_PROGRAMS=test_springreader test_vector \
			   test_linear_spring test_torsion_spring \
			   test_rk4 test_model test_residue \
			   test_sterics
TESTS=test_springreader test_vector \
	  test_linear_spring test_torsion_spring \
	  test_rk4 test_model test_residue \
	  test_sterics
BUILT_SOURCES=data/AA.c
SUFFIXES=.gperf

test_springreader_CFLAGS=$(OPENMP_CFLAGS)
test_springreader_SOURCES=t/springreader.c t/tap.c src/springreader.c \
						src/model.c src/residue.c src/sterics.c \
						data/AA.c src/linear_spring.c src/torsion_spring.c \
						src/vector.c

test_vector_CFLAGS=$(OPENMP_CFLAGS)
test_vector_SOURCES=t/vector.c src/vector.c t/tap.c

test_linear_spring_CFLAGS=$(OPENMP_CFLAGS)
test_linear_spring_SOURCES=t/linear_spring.c src/linear_spring.c \
						   src/sterics.c \
						   src/vector.c src/residue.c data/AA.c t/tap.c

test_torsion_spring_CFLAGS=$(OPENMP_CFLAGS)
test_torsion_spring_SOURCES=t/torsion_spring.c src/torsion_spring.c \
							src/sterics.c \
							src/vector.c src/residue.c data/AA.c t/tap.c

test_rk4_CFLAGS=$(OPENMP_CFLAGS)
test_rk4_SOURCES=t/rk4.c src/rk4.c src/model.c src/residue.c \
				 src/linear_spring.c src/torsion_spring.c \
				 src/springreader.c src/sterics.c \
				 src/vector.c data/AA.c t/tap.c

test_model_CFLAGS=$(OPENMP_CFLAGS)
test_model_SOURCES=t/model.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/sterics.c \
				   data/AA.c t/tap.c

test_residue_CFLAGS=$(OPENMP_CFLAGS)
test_residue_SOURCES=t/residue.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/springreader.c \
				   src/sterics.c data/AA.c t/tap.c

test_sterics_CFLAGS=$(OPENMP_CFLAGS)
test_sterics_SOURCES=t/sterics.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/springreader.c \
				   src/sterics.c data/AA.c t/tap.c

.gperf.c:
	$(GPERF) $< --output-file $@
