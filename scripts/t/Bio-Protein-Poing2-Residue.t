#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 8;

BEGIN {use_ok 'Bio::Protein::Poing2::Residue'}

my $res1 = Bio::Protein::Poing2::Residue->new(index => 1, type => 'A');
my $res2 = Bio::Protein::Poing2::Residue->new(index => 2, type => 'ALA');

is($res1->threeletter, 'ALA', 'One   -> Three letters');
is($res2->oneletter,   'A',   'Three -> One letters');
is($res2->index,       2,     'Index set');
is("$res1",            'ALA', 'Stringification works');
is_deeply($res1->atoms, [],   'Defaults with no atoms');

$res1->add_sidechain();
is($res1->atoms->[0]->name, 'ALA', 'Added sidechain atom');

my $gly = Bio::Protein::Poing2::Residue->new(index => 2, type => 'GLY');
$gly->add_sidechain();
is(@{$gly->atoms}, 0, 'No sidechain for glycine');
