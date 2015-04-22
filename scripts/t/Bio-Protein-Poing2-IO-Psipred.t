#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 19;
use File::Temp;

BEGIN {use_ok 'Bio::Protein::Poing2::IO::Psipred'}

my $ss = <<EOSS2;
# Confidences are dummy values
# File generated from hhblits MSA
   1 V C   1.000  0.000  0.000
   2 F C   1.000  0.000  0.000
   3 K E   0.000  0.000  1.000
   4 K E   0.000  0.000  1.000
   5 V E   0.000  0.000  1.000
   6 L E   0.000  0.000  1.000
   7 L E   0.000  0.000  1.000
   8 T E   0.000  0.000  1.000
EOSS2

open my $ss_in, q{<}, \$ss;
my $res1 = Bio::Protein::Poing2::IO::Psipred::read_psipred($ss_in);
close $ss_in;

is(keys %{$res1}, 8, 'Read 8 residues');
is($res1->{$_}->ss_state, 'C', "Residue $_ is C") for (1 .. 2);
is($res1->{$_}->ss_state, 'E', "Residue $_ is E") for (3 .. 8);

my $tmp = File::Temp->new;
print {$tmp} $ss;
$tmp->flush;

my $res2 = Bio::Protein::Poing2::IO::Psipred::read_psipred("$tmp");
is(keys %{$res2}, 8, 'Read 8 residues');
is($res2->{$_}->ss_state, 'C', "Residue $_ is C") for (1 .. 2);
is($res2->{$_}->ss_state, 'E', "Residue $_ is E") for (3 .. 8);
