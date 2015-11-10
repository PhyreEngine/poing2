package Bio::Protein::Poing2::LinearSpring;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use namespace::autoclean;
use Moose;

=head1 NAME

Bio::Protein::Poing2::LinearSpring - Connection between two atoms

=head1 SYNOPSIS

    #Spring between $atom1 and $atom2 with an equilibrium distance of 10A.
    my $spring = Bio::Protein::Poing2::LinearSpring->new(
        atom_1 => $atom1,
        atom_2 => $atom2,
        distance => 10,
    );
    print $spring->distance;

=head1 METHODS

=over

=item C<new(%args)> Instantiate a new class of type C<$type>.

Requires the arguments C<atom_1>, C<atom_2> and C<distance>.

=cut

has atom_1 => (is => 'ro', required => 1);
has atom_2 => (is => 'ro', required => 1);

#Inner and outer are used to compute handedness
has handedness => (is => 'ro', lazy => 1, builder => '_build_handedness', init_arg => undef);
has inner_atom => (is => 'rw', isa => 'Maybe[Bio::Protein::Poing2::Atom]', required => 0);
has outer_atom => (is => 'rw', isa => 'Maybe[Bio::Protein::Poing2::Atom]', required => 0);

has distance => (is => 'ro', required => 1);
has constant => (is => 'rw', default => undef);
has cutoff => (is => 'rw', default => undef);

sub string_repr {
    my ($self) = @_;
    return sprintf "% 4d % 4s % 4d % 4s %8.6f\n",# %8.6f %8.6f\n",
        $self->atom_1->residue->index, $self->atom_1->name,
        $self->atom_2->residue->index, $self->atom_2->name,
        $self->distance;
}

=item C<handedness()>

Calculate handedness (left of right) of this spring, using the C<inner> and
C<outer> atoms as the inner pair.

=cut

sub _build_handedness {
    my ($self) = @_;

    my $a = $self->atom_1;
    my $b = $self->atom_2;
    my $inner = $self->inner_atom;
    my $outer = $self->outer_atom;
    return undef unless $inner && $outer;

    my $ab = $b->coords     - $a->coords;
    my $ai = $inner->coords - $a->coords;
    my $ao = $outer->coords - $a->coords;
    my $cross = $ab x $ai;
    my $dot = $cross * $ao;
    return ($dot > 0) ? "RIGHT" : "LEFT";
}

=item C<TO_JSON()>

Serialise into a format usable by the JSON config file.

=cut

sub TO_JSON {
    my ($self) = @_;
    my $repr = {
        atoms => [$self->atom_1->index, $self->atom_2->index],
        distance => $self->distance,
    };
    $repr->{constant} = $self->constant if $self->constant;
    $repr->{cutoff}   = $self->cutoff   if $self->cutoff;
    $repr->{handedness} = {
        inner => $self->inner_atom->index,
        outer => $self->outer_atom->index,
        handedness => $self->handedness,
    } if $self->handedness;
    return $repr;
}

=back

=cut

__PACKAGE__->meta->make_immutable;
1;
