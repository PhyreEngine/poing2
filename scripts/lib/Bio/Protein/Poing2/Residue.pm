package Bio::Protein::Poing2::Residue;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data;

use Bio::Protein::Poing2::Class;
use if $^V gt v5.10.1, parent => 'Bio::Protein::Poing2::Class';
use if $^V le v5.10.1, base   => 'Bio::Protein::Poing2::Class';


#Allow overloading to string
use overload q{""} => 'threeletter';

=head1 NAME

Bio::Protein::Poing2::Residue - Class representing a residue

=head1 SYNOPSIS

    my $res1 = Bio::Protein::Poing2::Residue->new('ALA');
    my $res2 = Bio::Protein::Poing2::Residue->new('A');

    print $res1->oneletter;   #A
    print $res1->threeletter; #ALA
    print "$res1";            #ALA

=head1 METHODS

=over

=item C<new(type => $type)> Instantiate a new class of type C<$type>.

=cut

has type => (is => 'ro', required => 1);

=item C<atoms()>: Returns an arrayref of L<Bio::Protein::Poing2::Atom> objects.

=cut

has atoms => (is => 'ro', default => sub{[]});

=item C<ss_state()>: Returns a single-letter code for the
secondary structure state.

=cut

has ss_state => (is => 'ro');

=item C<oneletter()>: Returns the one-letter code for this residue.

=cut

sub oneletter {
    my ($self) = @_;
    my $aa = $self->type;
    $aa = $Bio::Protein::Poing2::Data::three2one{$aa} if length $aa == 3;
    return $aa;
}

=item C<threeletter()>: Returns the three-letter code for this residue.

=cut

sub threeletter {
    my ($self) = @_;
    my $aa = $self->type;
    $aa = $Bio::Protein::Poing2::Data::one2three{$aa} if length $aa == 1;
    return $aa;
}

=back

=cut
1;
