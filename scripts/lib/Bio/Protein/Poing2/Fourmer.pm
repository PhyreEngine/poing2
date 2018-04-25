package Bio::Protein::Poing2::Fourmer;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Data qw($PI);
use Math::Vector::Real;
use Moose::Role;

=head1 NAME

Bio::Protein::Poing2::Fourmer - Base class for torsion constraints on dihedrals.

=head1 SYNOPSIS

    my $fourmer = Bio::Protein::Poing2::Fourmer->new(
        atoms => [$atom_1, $atom_2, $atom_3, $atom_4],
    );
    print $fourmer->dihedral;

=cut

has atoms => (is => 'ro', required => 1);
has constant => (is => 'ro', isa => 'Maybe[Num]', default => undef);

requires 'dihedral';

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

1;
