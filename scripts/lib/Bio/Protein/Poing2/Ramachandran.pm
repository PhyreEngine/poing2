package Bio::Protein::Poing2::Ramachandran;
use strict;
use warnings;
use utf8;
use Moose;

use constant {
    PROLINE     => 'PROLINE',
    GLYCINE     => 'GLYCINE',
    ALANINE     => 'ALANINE',
    GENERAL     => 'GENERAL',
    PRE_PROLINE => 'PRE_PROLINE',
    ALPHA       => 'ALPHA',
    BETA        => 'BETA',
};

=head1 NAME

Bio::Protein::Poing2::Ramachandran - Represent a Ramachandran constraint

=head1 SYNOPSIS

    my $rama = Bio::Protein::Poing2::Ramachandran->new(
        residue => $residue,
        type    => Bio::Protein::Poing2->GENERAL,
    );
    print $rama->type;
    print $rama->string_repr;

=head1 SEE ALSO

L<Bio::Protein::Poing2::Ramachandran::List>

=cut

has residue => (
    is       => 'ro',
    required => 1,
    isa      => 'Bio::Protein::Poing2::Residue'
);
has type => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build_type'
);

use overload q{""} => \&string_repr;

sub _build_type {
    my ($self) = @_;

    #Base state off of secondary structure prediction if supplied, or residue
    #type if not.
    if($self->residue->ss_state && $self->residue->ss_state ne 'C'){
        if($self->residue->ss_state eq 'H'){
            return ALPHA;
        }elsif($self->residue->ss_state eq 'E'){
            return BETA;
        }
    }else{
        if($self->residue->threeletter eq 'PRO'){
            return PROLINE;
        }elsif($self->residue->threeletter eq 'GLY'){
            return GLYCINE;
        }elsif($self->residue->threeletter eq 'ALA'){
            return ALANINE;
        }else{
            return GENERAL;
        }
    }
}

sub string_repr {
    my ($self) = @_;
    return sprintf "% 4d %s\n", $self->residue->index, $self->type;
}

__PACKAGE__->meta->make_immutable;
1;
