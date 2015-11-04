package Bio::Protein::Poing2::Handedness;
use strict;
use warnings;
use utf8;
use Carp;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Handedness - Is a linear spring right or left handed?

=head1 SYNOPSIS

    my $handedness = Bio::Protein::Poing2::Handedness->new(
        atoms => [$a1, $a2, $a3, $a4],
    );
    print $handedness->string_repr;
    #10 13 17 20 RIGHT

=head1 ATTRIBUTES

=over

=item B<atoms> Arrayref of atoms.

=cut

has atoms => (is => 'rw', default => sub{[]});

sub handedness {
    my ($self) = @_;

    my $ab = $self->atoms->[3]->position - $self->atoms->[0]->position;
    my $ai = $self->atoms->[1]->position - $self->atoms->[0]->position;
    my $ao = $self->atoms->[2]->position - $self->atoms->[0]->position;
    my $cross = $ab x $ai;
    my $dot = $cross * $ao;
    return ($dot > 0) ? "RIGHT" : "LEFT";
}


sub string_repr {
    my ($self) = @_;
    my @residue_idx = map {$_->residue->index} @{$self->atoms};
    return join(q{ }, @residue_idx, $self->handedness), "\n";
}

__PACKAGE__->meta->make_immutable;
1;
