package Bio::Protein::Poing2::Atom;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data;
use Bio::Protein::Poing2;
use Moose;
use Math::Vector::Real;


=head1 NAME

Bio::Protein::Poing2::Atom - An atom

=head1 SYNOPSIS

    my $atom = Bio::Protein::Poing2::Atom->new(name => 'C');
    $atom->x(0.1);
    $atom->y(0.3);
    $atom->coords([0.1, 0.3, 0.5]);

=head1 METHODS

=over

=item C<new($name)>: Create an atom object. C<$name> can be anything, such as
C, CA, N, O, CB, etc.

=item C<name()>: Get the name of the atom.

=cut

has name   => (is => 'ro', required => 1);

=item C<index(): The ID of the atom.

=cut

has index => (is => 'rw', required => 1);


=item C<coords([$coords])>: Get the coordinates as an arrayref of (x, y, z)
values or set them to the values in C<$coords> if supplied.

=cut

has coords => (is => 'rw', default => sub { V(0, 0, 0) }
);

=item C<residue([$res]): Get/set the residue of the atom. This is optional,
because an atom does not necessarily belong to a residue.

=cut

has residue => (is => 'ro', default => undef);

=item C<x([$x])>: Get the x coordinate, or set the x coordinate to C<$x> if
supplied.

=cut

#Nobody is going to mistake ->x() for "a" x 10
##no critic (Subroutines::ProhibitBuiltinHomonyms)
sub x{
##use critic
    my ($self, $x) = @_;
    $self->coords->[0] = $x if defined $x;
    return $self->coords->[0];
}

=item C<y([$y])>: Get the y coordinate, or set the y coordinate to C<$y> if
supplied.

=cut

#Likewise, nobody will mistake ->y() for y///
##no critic (Subroutines::ProhibitBuiltinHomonyms)
sub y{
##use critic
    my ($self, $y) = @_;
    $self->coords->[1] = $y if defined $y;
    return $self->coords->[1];
}

=item C<z([$z])>: Get the z coordinate, or set the z coordinate to C<$z> if
supplied.

=cut

sub z{
    my ($self, $z) = @_;
    $self->coords->[2] = $z if defined $z;
    return $self->coords->[2];
}

=item C<string_repr()>: Get a string representation for the config file.

=cut

our $atom_rec = "ATOM  % 5d  %s %s A% 4d    %8.3f%8.3f%8.3f\n";
sub string_repr {
    my ($self) = @_;
    return sprintf $atom_rec,
        $self->index,
        sprintf("%-3s", $self->name),
        $self->residue->threeletter, $self->residue->index,
        $self->x, $self->y, $self->z;
}

=back

=cut

if(defined __PACKAGE__->meta){
    __PACKAGE__->meta->make_immutable;
}
1;
