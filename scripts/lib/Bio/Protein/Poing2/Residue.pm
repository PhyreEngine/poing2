package Bio::Protein::Poing2::Residue;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data;

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

=item C<new($type)> Instantiate a new class of type C<$type>.

=cut

sub new {
    my ($class, $type) = @_;

    croak "No residue type supplied" unless $type;
    my $aa = undef;
    if(length $type == 1){
        $aa = $Bio::Protein::Poing2::Data::one2three{$type};
    }elsif(length $type == 3){
        $aa = $type;
    }else{
        croak "Unknown amino acid code '$type': Must be one or three letters.";
    }

    return bless {
        type  => $aa,
        atoms => [],
    }, $class;
}

=item C<oneletter()>: Returns the one-letter code for this residue.

=cut

sub oneletter {
    my ($self) = @_;
    return $Bio::Protein::Poing2::Data::three2one{$self->threeletter};
}

=item C<threeletter()>: Returns the three-letter code for this residue.

=cut

sub threeletter {
    my ($self) = @_;
    return $self->{type};
}

=item C<atoms()>: Returns an arrayref of L<Bio::Protein::Poing2::Atom> objects.

=cut

sub atoms {
    my ($self) = @_;
    return $self->{atoms};
}

=back

=cut
1;
