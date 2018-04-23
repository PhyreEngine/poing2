# Poing2

Poing2 is a tool for combining multiple models of a protein. Constraints
describing pairwise residue distances and torsion angles are read from models
and combined with a force field constraining 3-atom angles, omega torsion
angles, and allowed regions of Ramachandran space.

## Setup

To compile poing2 from git, you will need [autoconf] and [automake]. First,
initialise the [cJSON] submodule:
```
git submodule init
git submodule update
```

Then, generate the `configure` script and `Makefile` with `autoreconf`. The
`-v` switch enables verbose mode, and the `-i` switch tells `autoreconf` to
install any missing files.
```
autoreconf -vi
```

To build the `poing2` executable, run the standard `./configure` and `make`
commands. Note that it is currently *not* recommended that you use OpenMP, so
disable it explicitly:
```
./configure --disable-openmp && make
```
Running `make` will require the [gperf] executable, which is used to generate a
perfect hash table for the atom types used by poing2.

Poing2 is packaged with scripts to generate JSON-formatted configuration files
from PDB-formatted models and FASTA alignments. These supporting scripts are
written in Perl, and require the following modules to be installed from CPAN:
[Moose][moose], [JSON::XS][json-xs], [List::BinarySearch][list-binarysearch],
and [Math::Vector::Real][math-vector-real].

Finally, Poing2 requires the `RAMA_DATA` environment variable to be set to the
directory containing files describing the allowed regions of Ramachandran
space. This is the `data` directory, located in the root of the project
alongside this `README.md` file, although you may move it wherever you please.
```
export RAMA_DATA=$PWD/data
```

You may also wish to set the `PERL5LIB` environment variable to point to the
`scripts/lib` directory so you can run the included scripts without the `-I`
switch. This document assumes that has not been done, and will include the `-I`
option in any sample commands.

## Running

Let's say you've been to the [Phyre2][phyre2] web server and got some models
for your query sequence. Your query sequence is stored in `query.fasta` and the
models are named `01.pdb`, `02.pdb`, etc. The alignments between the query
sequence and the PDB files are stored in correspondingly-named FASTA files:
`01.fasta`, `02.fasta`, etc.

**Note**: The *entire* query sequence must be present in each alignment file.
It is *not* good enough to trim the flanking regions.

**Note**: The residues in the PDB-like files containing each model *must* be
numbered consecutively, with the first *query* residue being residue 1. For
example, if the alignment of your model covers residues 10-50 of the query, the
residues in that PDB file must be numbered from 10-50. Models produced by
[Phyre2][phyre2] fit these criteria.

Generate a JSON configuration file using `scripts/bin/build_config.pl`:
```
perl -Iscripts/lib scritps/bin/build_config.pl \
    query.fasta \
    -t 01.fasta=01.pdb \
    -t 02.fasta=02.pdb \ # etc
    > config.json
```
    
The `build_config.pl` script takes lots of options, which can be seen by
running the script with the `--help` option.

The recommended parameters are:

+ **`--synth-time` 50**: Reduce the number of time-steps between synthesising
  atoms.  This will speed up the simulation greatly compared to the default
  value of 500 time steps.

+ **`--no-water`**: Disable the water bombardment simulation. Testing shows that
  the water bombardment tends not to show any benefits.

+ **`--no-shield-drag`**: Disable the shielded drag model, in which drag is not
  applied to atoms closer than steric radius plus the radius of gyration of a
  water molecule. The cruder drag model used when this is disabled results in
  models that are at least as accurate, and allows the simulation to converge
  much more quickly.

After generating the `config.json` file, run `poing2`:
```
./poing2 -s 100 config.json > model.pdb
```

The `-s 100` option tells poing2 to take a snapshot of the system every 100
time steps. Snapshots are written to standard output, and separated by
`MODEL`/`ENDMDL` records. Side-chains in Poing2 are modelled by a single
sphere, the centre of the sphere is recorded in the PDB file as an `ATOM`
record with the atom type corresponding to the residue type. That is, the
centre of the side-chain sphere for a `VAL` residue will have an atom type
of `VAL`.

[autoconf]: https://www.gnu.org/software/autoconf/autoconf.html
[automake]: https://www.gnu.org/software/automake/
[gperf]:https://www.gnu.org/software/gperf/
[cJSON]: https://github.com/DaveGamble/cJSON
[math-vector-real]: http://search.cpan.org/~salva/Math-Vector-Real-0.17/lib/Math/Vector/Real.pm
[list-binarysearch]: http://search.cpan.org/~davido/List-BinarySearch-0.25/lib/List/BinarySearch.pm
[json-xs]: http://search.cpan.org/~mlehmann/JSON-XS-3.02/XS.pm
[moose]: http://search.cpan.org/~ether/Moose-2.2010/lib/Moose.pm
[phyre2]: http://www.sbg.bio.ic.ac.uk/phyre2/
