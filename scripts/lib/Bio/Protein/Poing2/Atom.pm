package Bio::Protein::Poing2::Atom;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data;

=head1 NAME

Bio::Protein::Poing2::Atom - An atom

=head1 SYNOPSIS

    my $atom = Bio::Protein::Poing2::Atom->new('C');
    $atom->x(0.1);
    $atom->y(0.3);
    $atom->coords([0.1, 0.3, 0.5]);

=head1 METHODS

=over

=item C<new($name)>: Create an atom object. C<$name> can be anything, such as
C, CA, N, O, CB, etc.

=cut

sub new {
    my ($class, $name) = @_;
    croak "No name supplied" unless $name;

    return bless {
        name   => $name,
        coords => [0, 0, 0],
    }, $class;
}

=item C<name()>: Get the name of the atom.

=cut

sub name {
    my ($self) = @_;
    return $self->{name};
}

=item C<x([$x])>: Get the x coordinate, or set the x coordinate to C<$x> if
supplied.

=cut

sub x{
    my ($self, $x) = @_;
    $self->{coords}[0] = $x if defined $x;
    return $self->{coords}[0];
}

=item C<y([$y])>: Get the y coordinate, or set the y coordinate to C<$y> if
supplied.

=cut

sub y{
    my ($self, $y) = @_;
    $self->{coords}[1] = $y if defined $y;
    return $self->{coords}[1];
}

=item C<z([$z])>: Get the z coordinate, or set the z coordinate to C<$z> if
supplied.

=cut

sub z{
    my ($self, $z) = @_;
    $self->{coords}[2] = $z if defined $z;
    return $self->{coords}[2];
}

=item C<coords([$coords])>: Get the coordinates as an arrayref of (x, y, z)
values or set them to the values in C<$coords> if supplied.

=cut

sub coords {
    my ($self, $coords) = @_;
    $self->{coords} = $coords if $coords;
    return $self->{coords};
}

=back

=cut
1;
