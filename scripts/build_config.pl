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
use Bio::Protein::Poing2::Filter::Pair::MaxDistance;
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

=item B<--synth-time> I<T>

Wait I<T> time units between each atom synthesis Default: 500

=item B<--until> I<T>

Run for I<T> time units. Default: synth-time x (length + 5)

=item B<--no-sterics>

Do not steric forces. Default: Use sterics.

=item B<--no-shield-drag>

Use simplified drag force. Default: Do not simplify.

=item B<--no-water>

Do not bombard model. Default: Do bombardment.

=item B<--bb-only>

Only generate backbone atoms, not sidechains.

=item B<--fix-before> (Default: 50, disabled if B<--record-jitter> is false)

Requires this many atoms, from the most-recently synthesised atom backwards, to
remain free. The previous atoms may be fixed if their jitter becomes low enough.

=item B<--record-time> (Default: 10 × timestep, disabled if B<--record-jitter> is false)

Record the jitter at this time interval. Can be useful to avoid recording
millions of points when operating with a large C<fix_before> and and
C<synth_time> value.

=item B<--max-jitter> (Default: 0.01 Å)

Atoms with average jitter below this value are frozen if the B<--record-jitter>
option is set.

=item B<--record_jitter> (Default: false)

Should jitter be recorded and atoms near equilibrium be frozen? If true, the
B<--fix-before> B<record-time> and B<max-jitter> options are stored in the
configuration.

=item B<--max-distance> I<DISTANCE>

Keep springs with a distance less than or equal to I<DISTANCE> Angstroms.

=item B<--explicit-ss>

Add explicit springs and torsion constraints for predicted SS elements.

=item B<--position-from> I<FILE>

Read atom positions from the PDB file I<FILE>.

=item B<--add-hbonds>

Run stride on each template and add H bond springs and angles.

=item B<--hydrophobic-springs>

Add springs between each pair of hydrophobic residues.

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
    'synth-time=f',
    'until=f',
    'no-shield-drag',
    'no-water',
    'no-sterics',
    'bb-only',

    'fix-before=i',
    'record-time=f',
    'max-jitter=f',
    'record-jitter',

    'max-distance=f',
    'explicit-ss',
    'add-hbonds',
    'hydrophobic-springs',

    'position-from=s',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

my $query_file = shift || pod2usage('No sequence provided.');

my $query = Bio::Protein::Poing2::Query->new(
    sequence => $query_file,
    ss       => $options{ss},
    $options{'bb-only'}     ? (bb_only => 1)     : (),
    $options{'explicit-ss'} ? (explicit_ss => 1) : (),
    $options{'hydrophobic-springs'} ? (hydrophobic_springs => 1) : (),
    $options{'position-from'} ? (positions_file => $options{'position-from'}) : (),
);
my @templates = map {
    my ($aln, $model) = split q{=};
    Bio::Protein::Poing2::Template->new(
        alignment => $aln,
        model     => $model,
        query     => $query,
        add_hbonds => exists $options{'add-hbonds'},
)} @{$options{template}};

my @filters = (
    Bio::Protein::Poing2::Filter::Pair::SeqSep->new(min_sep => 3),
);
if($options{'max-distance'}){
    push @filters, Bio::Protein::Poing2::Filter::Pair::MaxDistance->new(
        distance => $options{'max-distance'},
    );
}

my %poing2_args = (
    query     => $query,
    templates => \@templates,
    verbose   => $options{verbose},
    $options{'synth-time'}  ? (synth_time  => $options{'synth-time'}) : (),
    $options{'until'}       ? (until       => $options{'until'}     ) : (),
    $options{'no-water'}    ? (use_water   => 0                     ) : (),
    $options{'no-sterics'}  ? (use_sterics => 0                     ) : (),
    $options{'shield_drag'} ? (shield_drag => 0                     ) : (),
    spring_filters => \@filters,
);
if($options{'record-jitter'}){
    $poing2_args{fix_before}    = $options{'fix-before'}  if $options{'fix-before'};
    $poing2_args{record_time}   = $options{'record-time'} if $options{'record-time'};
    $poing2_args{max_jitter}    = $options{'max-jitter'}  if $options{'max-jitter'};
    $poing2_args{record_jitter} = 1;
}

my $poing2 = Bio::Protein::Poing2->new(%poing2_args);

my $encoder = JSON::XS->new->convert_blessed->pretty;
print $encoder->encode($poing2);
