#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use autodie;
use Pod::Usage;
use Getopt::Long;
use Bio::Protein::Poing2;
use Bio::Protein::Poing2::Query;
use Bio::Protein::Poing2::Template;
use Bio::Protein::Poing2::Filter::Pair::SeqSep;
use JSON::XS;

=head1 NAME

build_config.pl - Build confing for poing2 run

=head1 USAGE

B<build_config.pl> B<query.fasta> [B<-t> I<aln.fasta>=I<model.pdb>] ...

=head1 OPTIONS

=over

=item B<-t>, B<--template> I<template.fasta>=I<template.pdb>

Add a template with alignment I<template.fasta> and model I<template.pdb>. This
option can be repeated as many times as necessary.

=item B<-s>, B<--ss> I<psipred.ss2>

Psipred SS2 file containing a secondary structure prediction.

=item B<-v>, B<--verbose>

Be a little more verbose about what the script is doing.

=item B<-h>, B<--help>

Display this help text.

=back

=cut

my %options = (
    verbose => 0,
);

Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
    'template|t=s@',
    'verbose|v',
    'ss|s=s',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

my $query_file = shift || pod2usage('No sequence provided.');

print STDERR "Parsing query\n" if $options{verbose};
my $query = Bio::Protein::Poing2::Query->new(
    sequence => $query_file,
    ss       => $options{ss},
);
my @templates = map {
    my ($aln, $model) = split q{=};
    print STDERR "Parsing $model\n" if $options{verbose};
    Bio::Protein::Poing2::Template->new(
        alignment => $aln,
        model     => $model,
        query     => $query,
)} @{$options{template}};

print STDERR "Combining templates and query\n" if $options{verbose};
my $poing2 = Bio::Protein::Poing2->new(
    query     => $query,
    templates => \@templates,
    spring_filters => [
        Bio::Protein::Poing2::Filter::Pair::SeqSep->new(min_sep => 3),
    ],
);
my $encoder = JSON::XS->new->convert_blessed->pretty;
print $encoder->encode($poing2);
