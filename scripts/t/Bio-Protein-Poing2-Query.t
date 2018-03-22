#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use Test::Exception;
use File::Temp;

BEGIN {use_ok 'Bio::Protein::Poing2::Query'}

my $alignment = <<EOS;
>T0644.fasta
DDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEP
STEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV
EOS

#my $model = do {local $/ = undef; <DATA>};

my $tmp_aln = File::Temp->new;
print {$tmp_aln} $alignment;
$tmp_aln->flush;

my $query = Bio::Protein::Poing2::Query->new(
    sequence => "$tmp_aln",
);

my @sequence = split //,
    'DDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPA'
    .'VLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV';
subtest "Sequence of query" => sub {
    plan(tests => 141);
    for(my $i=0; $i < @sequence; $i++){
        is($query->residues->{$i+1}->oneletter, $sequence[$i],
            "residue $i correct in query");
    }
}
