package Bio::Protein::Poing2::Filter::Pair::SeqSep;
use strict;
use warnings;
use Bio::Protein::Poing2;
use Moose;
extends 'Bio::Protein::Poing2::Filter::Pair';

=head1 NAME

Bio::Protein::Poing2::Filter::Pair::SeqSep - Filter by sequence separation

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Pair::SeqSep->new(
        max_sep => 10,
        min_sep => 5,
    );

    my $r1 = Bio::Protein::Poing2::Residue->new(index => 1,  type => 'A');
    my $r2 = Bio::Protein::Poing2::Residue->new(index => 20, type => 'G');

    my $a1 = Bio::Protein::Poing2::Atom->new(name => 'CA', residue => $r1);
    my $a2 = Bio::Protein::Poing2::Atom->new(name => 'CA', residue => $r2);

    my $pair = Bio::Protein::Poing2::LinearSpring->new(
        atom_1 => $a1,
        atom_2 => $a2,
    );
    my $pairs = $filt->filter([$pair]);

=head1 DESCRIPTION

Filter the list of pairs (L<Bio::Protein::Poing2::LinearSpring> objects) such
that the sequence separation between each atom is less than or equal to
C<max_sep> residues, and greater than or equal to C<min_sep> residues.

Set C<max_sep> or C<min_sep> to C<undef> to allow an arbitrary separation.

=cut

has max_sep => (is => 'rw', required => 0);
has min_sep => (is => 'rw', required => 0);

sub filter {
    my ($self, $pairs) = @_;

    my @new_pairs = grep {
        my $r1 = $_->atom_1->residue->index;
        my $r2 = $_->atom_2->residue->index;
        my $sep = abs($r1 - $r2);
        (!defined($self->{min_sep}) || $sep >= $self->{min_sep})
        &&
        (!defined($self->{max_sep}) || $sep <= $self->{max_sep})
    } @{$pairs};
    return \@new_pairs;
}

__PACKAGE__->meta->make_immutable;

1;
