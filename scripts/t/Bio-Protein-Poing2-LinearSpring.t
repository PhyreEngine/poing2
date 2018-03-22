#!/usr/bin/env perl
use strict;
use warnings;
use Math::Vector::Real;
use Test::More tests => 2;
use Bio::Protein::Poing2::Atom;

BEGIN {use_ok 'Bio::Protein::Poing2::LinearSpring'}

my $spring = Bio::Protein::Poing2::LinearSpring->new(
    atom_1 => Bio::Protein::Poing2::Atom->new(
        index => 1,
        name => 'C',
        coords => V(0, 0, 0),
    ),
    atom_2 => Bio::Protein::Poing2::Atom->new(
        index => 2,
        name => 'C',
        coords => V(0, 1, 0),
    ),
    distance => 1.0,
);
cmp_ok($spring->distance, '==', 1.0, 'Distance correct');
