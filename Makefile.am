bin_PROGRAMS=poing2

poing2_deps=src/rk4.c src/leapfrog.c src/model.c src/residue.c \
			   src/linear_spring.c src/torsion_spring.c src/springreader.c \
			   src/bond_angle.c \
			   src/vector.c src/sterics.c data/atoms.c data/AA.c \
			   src/rama.c src/cJSON/cJSON.c src/rattle.c \
			   src/record.c
poing2_CFLAGS=$(OPENMP_CFLAGS)
poing2_SOURCES=src/poing.c $(poing2_deps)



check_PROGRAMS=test_springreader test_vector \
			   test_linear_spring test_torsion_spring \
			   test_model \
			   test_sterics test_bond_angle \
			   test_record
TESTS=test_springreader test_vector \
	  test_linear_spring test_torsion_spring \
	  test_model test_residue \
	  test_sterics test_bond_angle \
	  test_record

CLEANFILES=data/AA.c data/AA.h data/atoms.c data/atoms.h

if HAVE_CLOCK_GETTIME_AM
poing2_deps += src/profile.c
check_PROGRAMS += test_profile
TESTS += test_profile
test_profile_CFLAGS=$(OPENMP_CFLAGS)
test_profile_SOURCES=t/profile.c t/tap.c $(poing2_deps)
endif

test_record_CFLAGS=$(OPENMP_CFLAGS)
test_record_SOURCES=t/record.c t/tap.c $(poing2_deps)

test_springreader_CFLAGS=$(OPENMP_CFLAGS)
test_springreader_SOURCES=t/springreader.c t/tap.c $(poing2_deps)

test_vector_CFLAGS=$(OPENMP_CFLAGS)
test_vector_SOURCES=t/vector.c t/tap.c $(poing2_deps)

test_linear_spring_CFLAGS=$(OPENMP_CFLAGS)
test_linear_spring_SOURCES=t/linear_spring.c t/tap.c $(poing2_deps)

test_torsion_spring_CFLAGS=$(OPENMP_CFLAGS)
test_torsion_spring_SOURCES=t/torsion_spring.c t/tap.c $(poing2_deps)

test_bond_angle_CFLAGS=$(OPENMP_CFLAGS)
test_bond_angle_SOURCES=t/bond_angle.c t/tap.c $(poing2_deps)

test_model_CFLAGS=$(OPENMP_CFLAGS)
test_model_SOURCES=t/model.c t/tap.c $(poing2_deps)

test_sterics_CFLAGS=$(OPENMP_CFLAGS)
test_sterics_SOURCES=t/sterics.c t/tap.c $(poing2_deps)

data/atoms.c: data/atoms.gperf
	gperf $< --output-file $@
	sed -i 's/{""}/{"", 0, 0, 0, 0}/g' "$@"

data/AA.c: data/AA.gperf
	gperf $< --output-file $@
	sed -i 's/{""}/{"", "", 0, 0, 0, 0, 0}/g' "$@"
