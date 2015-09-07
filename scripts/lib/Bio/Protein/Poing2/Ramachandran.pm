package Bio::Protein::Poing2::Ramachandran;
use strict;
use warnings;
use utf8;
use Carp;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Ramachandran - Represent a Ramachandran constraint

=head1 SYNOPSIS

    my $rama = Bio::Protein::Poing2::Ramachandran->new(
        residue => $residue_index,
        type    => 'GENERAL',
    );
    print $rama->string_repr;

=cut

has residue => (is => 'ro', required => 1, isa => 'Int');
has type    => (is => 'ro', required => 1, isa => 'Str');

sub string_repr {
    my ($self) = @_;
    return sprintf "% 3d %s\n", $self->residue, $self->type;
}

__PACKAGE__->meta->make_immutable;
1;
