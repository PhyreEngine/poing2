#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 5;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;

BEGIN {use_ok 'Bio::Protein::Poing2::Filter::Atom::Backbone'}

my $filter = Bio::Protein::Poing2::Filter::Atom::Backbone->new;

my %residues = (
    1 => Bio::Protein::Poing2::Residue->new(
        index => 1,
        type => 'ALA',
        atoms => [
            Bio::Protein::Poing2::Atom->new(name => 'C', index => 1),
            Bio::Protein::Poing2::Atom->new(name => 'K', index => 2),
        ],
    ),
    2 => Bio::Protein::Poing2::Residue->new(
        index => 2,
        type => 'GLY',
        atoms => [
            Bio::Protein::Poing2::Atom->new(name => 'CA', index => 1),
            Bio::Protein::Poing2::Atom->new(name => 'U',  index => 2),
        ],
    ),
);

my $filt = $filter->filter(\%residues);
is($filt->{1}->atoms->[0]->name, 'C', 'Atom 1 is C');
is(@{$filt->{1}->atoms},           1, '1 atom in residue 1');

is($filt->{2}->atoms->[0]->name, 'CA', 'Atom 2 is CA');
is(@{$filt->{2}->atoms},            1, '1 atom in residue 2');
