#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use Bio::Protein::Poing2::Atom;

BEGIN {use_ok 'Bio::Protein::Poing2::LinearSpring'}

my $spring = Bio::Protein::Poing2::LinearSpring->new(
    atom_1 => Bio::Protein::Poing2::Atom->new(
        name => 'C',
        coords => Bio::Protein::Poing2::Vector->new(coords => [0, 0, 0]),
    ),
    atom_2 => Bio::Protein::Poing2::Atom->new(
        name => 'C',
        coords => Bio::Protein::Poing2::Vector->new(coords => [0, 1, 0]),
    ),
);
cmp_ok($spring->distance, '==', 1.0, 'Distance correct');
