#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 6;

use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::LinearSpring;
BEGIN {use_ok 'Bio::Protein::Poing2::Filter::Pair::SeqSep'}

my $filt = Bio::Protein::Poing2::Filter::Pair::SeqSep->new(
    max_sep => 10,
    min_sep => 5,
);

my $r1 = Bio::Protein::Poing2::Residue->new(index => 1,  type => 'A');
my $r2 = Bio::Protein::Poing2::Residue->new(index => 20, type => 'G');

my $a1 = Bio::Protein::Poing2::Atom->new(name => 'CA', residue => $r1);
my $a2 = Bio::Protein::Poing2::Atom->new(name => 'CA', residue => $r2);

my $pair = Bio::Protein::Poing2::LinearSpring->new(
    atom_1 => $a1,
    atom_2 => $a2,
    distance => 1.0,
);

is(@{$filt->filter([$pair])}, 0, 'Removed all springs');
$filt->max_sep(30);
is(@{$filt->filter([$pair])}, 1, 'Kept spring');
$filt->max_sep(undef);
is(@{$filt->filter([$pair])}, 1, 'Kept spring');
$filt->min_sep(25);
is(@{$filt->filter([$pair])}, 0, 'Removed all springs');
$filt->min_sep(undef);
is(@{$filt->filter([$pair])}, 1, 'Kept spring');
