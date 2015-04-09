package Bio::Protein::Poing2::Fourmer;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Vector;
use Bio::Protein::Poing2::Class;
use Bio::Protein::Poing2::Data qw($PI);
use base qw(Bio::Protein::Poing2::Class);

=head1 NAME

Bio::Protein::Poing2::Fourmer - Represent four atoms, so we can get dihedral angles

=head1 SYNOPSIS

    my $fourmer = Bio::Protein::Poing2::Fourmer->new(
        atoms => [$atom_1, $atom_2, $atom_3, $atom_4],
    );
    print $fourmer->dihedral;

=cut

has atoms => (is => 'ro', required => 1);

sub dihedral {
    my ($self) = @_;

    my $b1 = $self->atoms->[1]->coords - $self->atoms->[0]->coords;
    my $b2 = $self->atoms->[2]->coords - $self->atoms->[1]->coords;
    my $b3 = $self->atoms->[3]->coords - $self->atoms->[2]->coords;

    my $angle = atan2(
        ($b1 x $b2) x ($b2 x $b3) . $b2 / $b2->mag,
        ($b1 x $b2) . ($b2 x $b3)
    ) * 180 / $PI;
    return $angle;
}

1;
