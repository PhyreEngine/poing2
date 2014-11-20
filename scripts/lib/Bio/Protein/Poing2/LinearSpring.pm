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
has distance => (is => 'ro', required => 1);

=back

=cut
1;
