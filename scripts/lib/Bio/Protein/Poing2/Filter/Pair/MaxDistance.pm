package Bio::Protein::Poing2::Filter::Pair::MaxDistance;
use strict;
use warnings;
use Moose;
extends 'Bio::Protein::Poing2::Filter::Pair';

=head1 NAME

Bio::Protein::Poing2::Filter::Pair::MaxDistance - Discard sprints with equilibrium distance over a threshold.

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Pair::MaxDistance->new(
        distance => 10 #Angstroms
    );

=head1 DESCRIPTION

Keep springs that have an equilibrium distance below a certain value.

=cut

has distance => (is => 'ro', isa => 'Num', required => 1);

sub filter {
    my ($self, $pairs) = @_;

    my @new_pairs = ();
    for(@{$pairs}){
        if($_->distance <= $self->distance){
            push @new_pairs, $_;
        }
    }
    return \@new_pairs;
}

__PACKAGE__->meta->make_immutable;

1;
