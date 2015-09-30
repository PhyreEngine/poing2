package Bio::Protein::Poing2::Filter::Pair::ConsecutiveOnly;
use strict;
use warnings;
use Bio::Protein::Poing2;
use Bio::Protein::Poing2::Data qw(%fine_bb_links);
use Moose;

=head1 NAME

Bio::Protein::Poing2::Filter::Pair::ConsecutiveOnly - Keep pairs between consecutive atoms only.

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Pair::ConsecutiveOnly->new();

=head1 DESCRIPTION

Keep springs between consecutive atoms only. That is, only the following
springs will be kept:

         O
         ║
     N   C
    ╱ ╲ ╱ ╲
   C  CA   N
       │
       R

=cut

sub filter {
    my ($self, $pairs) = @_;

    my @new_pairs = ();
    for(@{$pairs}){
        my $a1 = $_->atom_1->name;
        my $a2 = $_->atom_2->name;
        my $r1 = $a1->residue->index;
        my $r2 = $a2->residue->index;

        #Link can be either way round; try both
        if(!exists $fine_bb_links{$a1}){
            my ($tmp_a, $tmp_r) = ($a1, $r1);
            $a1 = $a2;
            $a2 = $tmp_a;
            $r1 = $r2;
            $r2 = $tmp_r;
        }
        next unless exists $fine_bb_links{$a1};

        #Keep sidechain atoms
        if($a1 eq 'CA' && $Bio::Protein::Poing2::Data::CA_SC_len{$a2}){
            push @new_pairs, $_;
            next;
        }

        #Keep links
        for my $link(@{$fine_bb_links{$a1}}){
            next if $link->{atom} ne $a2;
            next if $r1 + $link->{increment} != $r2;
            push @new_pairs, $_;
        }
    }
    return \@new_pairs;
}

__PACKAGE__->meta->make_immutable;

1;
