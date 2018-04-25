package Bio::Protein::Poing2::Fourmer::Calculated;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data qw($PI);
use Moose;

=head1 NAME

Bio::Protein::Poing2::Fourmer::Calculated - Dihedral angle calculated from four atoms

=head1 SYNOPSIS

    my $fourmer = Bio::Protein::Poing2::Fourmer->new(
        atoms => [$atom_1, $atom_2, $atom_3, $atom_4],
        constant => 0.1, # Optional
    );
    print $fourmer->dihedral;

=head1 DESCRIPTION

This class can be used to represent a torsion constraint on the dihedral angle
between four atoms. The dihedral angle is calculated from the atoms, so this
class should be used when extracting torsion constraints from a template.

=cut

with 'Bio::Protein::Poing2::Fourmer';

sub dihedral {
    my ($self) = @_;

    my $b1 = $self->atoms->[1]->coords - $self->atoms->[0]->coords;
    my $b2 = $self->atoms->[2]->coords - $self->atoms->[1]->coords;
    my $b3 = $self->atoms->[3]->coords - $self->atoms->[2]->coords;

    my $angle = atan2(
        ($b1 x $b2) x ($b2 x $b3) * $b2 / abs($b2),
        ($b1 x $b2) * ($b2 x $b3)
    ) * 180 / $PI;
    $angle += 360 if $angle < 0;
    return $angle;
}

__PACKAGE__->meta->make_immutable;

1;
