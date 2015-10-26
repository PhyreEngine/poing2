#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 14;
use File::Temp;

BEGIN {use_ok 'Bio::Protein::Poing2::IO::FastaAlignment'}

my $fasta = <<EOFASTA;
>Query
AAA
GGG
>Template
GGG
-AA
EOFASTA

my $res;

#Test with glob
open my $fas_in, q{<}, \$fasta;
$res = Bio::Protein::Poing2::IO::FastaAlignment::read_fasta($fas_in);
close $fas_in;

is(keys %{$res}, 6, 'Read 6 residues successfully');
is($res->{1}->{from}->oneletter, 'A', "Residue q1");
is($res->{2}->{from}->oneletter, 'A', "Residue q2");
is($res->{3}->{from}->oneletter, 'A', "Residue q3");
is($res->{4}->{from}->oneletter, 'G', "Residue q4");
is($res->{5}->{from}->oneletter, 'G', "Residue q5");
is($res->{6}->{from}->oneletter, 'G', "Residue q6");

is($res->{1}->{to}->oneletter, 'G', "Residue t1");
is($res->{2}->{to}->oneletter, 'G', "Residue t2");
is($res->{3}->{to}->oneletter, 'G', "Residue t3");
is($res->{4}->{to},          undef, "Residue t4");
is($res->{5}->{to}->oneletter, 'A', "Residue t5");
is($res->{6}->{to}->oneletter, 'A', "Residue t6");
