package Bio::Protein::Poing2::Fourmer::Fixed;
use strict;
use warnings;
use utf8;
use Carp;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Fourmer::Fixed - A fixed torsion angle around four atoms

=head1 SYNOPSIS

    my $fourmer = Bio::Protein::Poing2::Fourmer::Fixed->new(
        atoms => [$atom_1, $atom_2, $atom_3, $atom_4],
        dihedral => 180,
        constant => 0.1, # Optional
    );
    print $fourmer->dihedral;

=head1 DESCRIPTION

This class is used to represent a constraint on the dihedral angle of a set of
four atoms. The dihedral angle is a fixed value, specified in degrees. An
optional force constant may also be supplied.

=cut

has dihedral => (is => 'ro', required => 1, isa => 'Num');
with 'Bio::Protein::Poing2::Fourmer';

__PACKAGE__->meta->make_immutable;

1;
