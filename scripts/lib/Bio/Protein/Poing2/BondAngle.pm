package Bio::Protein::Poing2::BondAngle;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use Math::Trig qw(acos pi);
use Moose;

=head1 NAME

Bio::Protein::Poing2::BondAngle - Angle between three atoms

=head1 SYNOPSIS

    #Angle with an equilibrium angle of 45 degrees
    my $spring = Bio::Protein::Poing2::BondAngle->new(
        atoms => [$atom1, $atom2, $atom3],
        angle => 45,
    );
    print $spring->angle;

=head1 METHODS

=over

=item C<new(%args)> Instantiate a new object.

Requires the argument C<atoms>.

=cut

has atoms => (is => 'ro', isa => 'ArrayRef', required => 1);
has constant => (is => 'rw', default => 0.1);
has cutoff => (is => 'rw', default => undef);
has angle => (is => 'ro', lazy => 1, builder => '_calculate_angle');

sub _calculate_angle {
    my ($self) = @_;
    my ($rij) = $self->atoms->[0]->coords - $self->atoms->[1]->coords;
    my ($rkj) = $self->atoms->[2]->coords - $self->atoms->[1]->coords;

    return acos( ($rij * $rkj) / (abs($rij) * abs($rkj)) ) * 180 / pi;
}

sub string_repr {
    my ($self) = @_;
    return sprintf "% 4d % 4s % 4d % 4s % 4d % 4s %8.6f %8.6f\n",
        $self->atoms->[0]->residue->index, $self->atoms->[0]->name,
        $self->atoms->[1]->residue->index, $self->atoms->[1]->name,
        $self->atoms->[2]->residue->index, $self->atoms->[2]->name,
        $self->angle, $self->constant;
}


=back

=cut

__PACKAGE__->meta->make_immutable;
1;
