bin_PROGRAMS=poing2

poing2_deps=src/rk4.c src/leapfrog.c src/model.c src/residue.c \
			   src/linear_spring.c src/torsion_spring.c src/springreader.c \
			   src/bond_angle.c \
			   src/vector.c src/sterics.c data/AA.c data/atoms.c \
			   src/rama.c src/cJSON/cJSON.c
poing2_CFLAGS=$(OPENMP_CFLAGS)
poing2_SOURCES=src/poing.c $(poing2_deps)


check_PROGRAMS=test_springreader test_vector \
			   test_linear_spring test_torsion_spring \
			   test_rk4 test_model test_residue \
			   test_sterics test_bond_angle
TESTS=test_springreader test_vector \
	  test_linear_spring test_torsion_spring \
	  test_model test_residue \
	  test_sterics test_bond_angle
BUILT_SOURCES=data/AA.c data/atoms.c
SUFFIXES=.gperf

test_springreader_CFLAGS=$(OPENMP_CFLAGS)
test_springreader_SOURCES=t/springreader.c t/tap.c $(poing2_deps)

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

test_bond_angle_CFLAGS=$(OPENMP_CFLAGS)
test_bond_angle_SOURCES=t/bond_angle.c src/bond_angle.c \
							src/sterics.c \
							src/vector.c src/residue.c data/AA.c t/tap.c

test_rk4_CFLAGS=$(OPENMP_CFLAGS)
test_rk4_SOURCES=t/rk4.c src/rk4.c src/model.c src/residue.c \
				 src/linear_spring.c src/torsion_spring.c src/bond_angle.c \
				 src/springreader.c src/sterics.c \
				 src/vector.c data/AA.c t/tap.c data/atoms.c

test_model_CFLAGS=$(OPENMP_CFLAGS)
test_model_SOURCES=t/model.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/sterics.c \
				   data/AA.c t/tap.c src/bond_angle.c

test_residue_CFLAGS=$(OPENMP_CFLAGS)
test_residue_SOURCES=t/residue.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/springreader.c \
				   src/sterics.c data/AA.c t/tap.c data/atoms.c src/bond_angle.c

test_sterics_CFLAGS=$(OPENMP_CFLAGS)
test_sterics_SOURCES=t/sterics.c src/model.c src/residue.c src/vector.c \
				   src/linear_spring.c src/torsion_spring.c src/springreader.c \
				   src/sterics.c data/AA.c t/tap.c data/atoms.c src/bond_angle.c \
				   src/rama.c

.gperf.c:
	$(GPERF) $< --output-file $@
