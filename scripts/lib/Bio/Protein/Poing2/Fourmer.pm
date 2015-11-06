package Bio::Protein::Poing2::Fourmer;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Data qw($PI);
use Math::Vector::Real;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Fourmer - Represent four atoms, so we can get dihedral angles

=head1 SYNOPSIS

    my $fourmer = Bio::Protein::Poing2::Fourmer->new(
        atoms => [$atom_1, $atom_2, $atom_3, $atom_4],
    );
    print $fourmer->dihedral;

=cut

has atoms => (is => 'ro', required => 1);

has constant => (is => 'ro', isa => 'Maybe[Num]', default => undef);

sub dihedral {
    my ($self) = @_;

    my $b1 = $self->atoms->[1]->coords - $self->atoms->[0]->coords;
    my $b2 = $self->atoms->[2]->coords - $self->atoms->[1]->coords;
    my $b3 = $self->atoms->[3]->coords - $self->atoms->[2]->coords;

    my $angle = atan2(
        ($b1 x $b2) x ($b2 x $b3) * $b2 / abs($b2),
        ($b1 x $b2) * ($b2 x $b3)
    ) * 180 / $PI;
    return $angle;
}

sub string_repr {
    my ($self) = @_;
    my $line = "% 4d %s % 4d %s % 4d %s % 4d %s %8.3f %8.3f\n";# %8.3f\n";
    return sprintf $line,
        $self->atoms->[0]->residue->index, $self->atoms->[0]->name,
        $self->atoms->[1]->residue->index, $self->atoms->[1]->name,
        $self->atoms->[2]->residue->index, $self->atoms->[2]->name,
        $self->atoms->[3]->residue->index, $self->atoms->[3]->name,
        $self->dihedral, $self->constant;
}

sub TO_JSON {
    my ($self) = @_;
    my $repr = {
        atoms => [map {$_->index} @{$self->atoms}],
        angle => $self->dihedral,
    };
    $repr->{constant} = $self->constant if $self->constant;
    return $repr;
}

__PACKAGE__->meta->make_immutable;

1;
