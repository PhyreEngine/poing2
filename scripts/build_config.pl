#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use Bio::Protein::Poing2;
use Bio::Protein::Poing2::IO::Fasta;
use Bio::Protein::Poing2::IO::PDB;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Filter::Atom::Backbone;
use Bio::Protein::Poing2::Filter::Pair::SeqSep;
use Bio::Protein::Poing2::Filter::Pair::PruneLong;

=head1 NAME

build_config.pl - Build confing for poing2 run

=head1 USAGE

B<build_config.pl> B<query.fasta> B<template1.pdb> [B<template2.pdb>] ...

=head1 ARGUMENTS

The first argument must be a FASTA file containing the query sequence. The
remaining arguments are assumed to be PDB files containing templates.

=head1 OPTIONS

=over

=item B<-h>, B<--help>

Display this help text.

=item B<--bb-only>

Only use backbone atoms.

=item B<-v>, B<--verbose>

Give a bit more information (to standard error) about what is happening.

=back

=cut

my %options = (
    verbose => 0,
);

Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
    'coarse',
    'verbose|v',
    'min-seq-sep=i',
    'max-seq-sep=i',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

my $query     = shift || pod2usage('No sequence provided.');
my @templates = @ARGV;

my @atom_filters = ();
my @pair_filters = ();
push @pair_filters, Bio::Protein::Poing2::Filter::Pair::SeqSep->new(
    min_sep => $options{'min-seq-sep'},
    max_sep => $options{'max-seq-sep'},
);

push @atom_filters, Bio::Protein::Poing2::Filter::Atom::Backbone->new();
push @pair_filters, Bio::Protein::Poing2::Filter::Pair::PruneLong->new();


my $poing2 = Bio::Protein::Poing2->new(
    query => $query,
    pdbs  => \@templates,
    verbose => $options{verbose},

    atom_filters => \@atom_filters,
    pair_filters => \@pair_filters,
);

if($options{coarse}){
    $poing2->init_coarse_bb;
    $poing2->init_coarse_sc;
}else{
    $poing2->init_fine_bb;
    $poing2->init_fine_bb_angles;
    $poing2->init_coarse_sc;
}
$poing2->renumber_atoms;
$poing2->set_positions;

print $poing2->string_repr;
