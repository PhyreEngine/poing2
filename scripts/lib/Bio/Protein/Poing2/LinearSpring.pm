package Bio::Protein::Poing2::LinearSpring;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Vector;
use Bio::Protein::Poing2::Class;
use base qw(Bio::Protein::Poing2::Class);

=head1 NAME

Bio::Protein::Poing2::LinearSpring - Connection between two atoms

=head1 SYNOPSIS

    my $spring = Bio::Protein::Poing2::LinearSpring->new(
        atom_1 => $atom1,
        atom_2 => $atom2,
    );
    print $spring->distance;

=head1 METHODS

=over

=item C<new(%args)> Instantiate a new class of type C<$type>.

Requires the arguments 'atom_1' and 'atom_2';

=cut

has atom_1 => (is => 'ro', required => 1);
has atom_2 => (is => 'ro', required => 1);

sub distance {
    my ($self) = @_;
    return ($self->{atom_1}->coords - $self->{atom_2}->coords)->mag;
}

=back

=cut
1;
