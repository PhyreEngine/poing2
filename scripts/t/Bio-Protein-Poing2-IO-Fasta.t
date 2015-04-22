#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 15;
use File::Temp;

BEGIN {use_ok 'Bio::Protein::Poing2::IO::Fasta'}

my $fasta = <<EOFASTA;
>Query
AAAGGG
EOFASTA

my $res;

#Test with glob
open my $fas_in, q{<}, \$fasta;
$res = Bio::Protein::Poing2::IO::Fasta::read_fasta($fas_in);
close $fas_in;

is(keys %{$res}, 6, 'Read 6 residues successfully');
is($res->{$_}->oneletter, 'A', "Residue $_ (glob)") for (1 .. 3);
is($res->{$_}->oneletter, 'G', "Residue $_ (glob)") for (4 .. 6);

my $tmpfile = File::Temp->new;
print {$tmpfile} $fasta;
$tmpfile->flush;

$res = Bio::Protein::Poing2::IO::Fasta::read_fasta("$tmpfile");
is(keys %{$res}, 6, 'Read 6 residues successfully');
is($res->{$_}->oneletter, 'A', "Residue $_ (file)") for (1 .. 3);
is($res->{$_}->oneletter, 'G', "Residue $_ (file)") for (4 .. 6);
