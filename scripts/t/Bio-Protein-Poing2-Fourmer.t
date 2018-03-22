#!/usr/bin/env perl
use strict;
use warnings;
use Math::Vector::Real;
use Test::More tests => 4;
use Bio::Protein::Poing2::Atom;

BEGIN {use_ok 'Bio::Protein::Poing2::Fourmer'}

my $atoms = [
    Bio::Protein::Poing2::Atom->new(
        index => 0,
        name => 'CA',
        coords => V(-1, 0, 0),
    ),
    Bio::Protein::Poing2::Atom->new(
        index => 1,
        name => 'CA',
        coords => V(0, 0, 0),
    ),
    Bio::Protein::Poing2::Atom->new(
        index => 2,
        name => 'CA',
        coords => V(0, 0, 1),
    ),
    Bio::Protein::Poing2::Atom->new(
        index => 3,
        name => 'CA',
        coords => V(-1, 0, 1),
    ),
];
my $fourmer = Bio::Protein::Poing2::Fourmer->new(
    atoms => $atoms,
);
is($fourmer->dihedral, 0, 'Dihedral angle is zero');

$atoms->[3] = Bio::Protein::Poing2::Atom->new(
    index => 3,
    name => 'CA',
    coords => V(1, 0, 1),
);
is($fourmer->dihedral, 180, 'Dihedral angle is 180 degrees');

$atoms->[3] = Bio::Protein::Poing2::Atom->new(
    index => 3,
    name => 'CA',
    coords => V(0, 1, 1),
);
is($fourmer->dihedral, 270, 'Dihedral angle is 270 degrees');
