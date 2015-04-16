package Bio::Protein::Poing2::Filter::Pair::PruneLong;
use strict;
use warnings;
use Bio::Protein::Poing2;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Filter::Pair::PruneLong - Keep only CA-CA springs for springs between different residues

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Pair::PruneLong->new();

=cut

sub filter {
    my ($self, $pairs) = @_;

    my @new_pairs = ();
    for(@{$pairs}){
        my $a1 = $_->atom_1;
        my $a2 = $_->atom_2;
        my $r1 = $a1->residue->index;
        my $r2 = $a2->residue->index;
        if($r1 == $r2 || ($a1->name eq 'CA' && $a2->name eq 'CA')){
            push @new_pairs, $_;
        }
    }
    return \@new_pairs;
}

__PACKAGE__->meta->make_immutable;

1;
