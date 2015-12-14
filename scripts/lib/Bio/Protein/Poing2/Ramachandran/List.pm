package Bio::Protein::Poing2::Ramachandran::List;
use strict;
use warnings;
use utf8;
use Moose;
use Bio::Protein::Poing2::Ramachandran;
use List::Util qw(min max);

=head1 NAME

Bio::Protein::Poing2::Ramachandran::List - Build a list of Ramachandran constraints

=head1 SYNOPSIS

    my $rama = Bio::Protein::Poing2::Ramachandran::List->new(
        residues => $residues,
    );
    print $rama->type;
    print $rama->string_repr;

=head1 SEE ALSO

L<Bio::Protein::Poing2::Ramachandran::List>

=cut

has residues => (
    is       => 'ro',
    required => 1,
    isa      => 'HashRef[Bio::Protein::Poing2::Residue]'
);

has ramachandrans => (
    is      => 'ro',
    lazy    => 1,
    builder => '_build_rama',
);

use overload q{""} => \&string_repr;

sub _build_rama {
    my ($self) = @_;

    my @rama = ();

    my $min = min keys %{$self->residues};
    my $max = max keys %{$self->residues};

    for my $index (sort {$a <=> $b} keys %{$self->residues}){
        #Skip the first and last residues, as we can't determine a phi/psi
        #angle for those
        next if $index == $min || $index == $max;

        my $resi      = $self->residues->{$index};
        my $next_resi = $self->residues->{$index + 1};

        my $rama = undef;

        #If the next residue is a proline, mark the current residue as PRE_PROLINE
        if($next_resi && $next_resi->threeletter eq 'PRO'){
            $rama = Bio::Protein::Poing2::Ramachandran->new(
                residue => $resi,
                type    => Bio::Protein::Poing2::Ramachandran->PRE_PROLINE,
            );
        }else{
            #Let the residue type determine the Ramachandran type by default
            $rama = Bio::Protein::Poing2::Ramachandran->new(
                residue => $resi,
            );
        }
        push @rama, $rama;
    }
    return \@rama;
}

sub string_repr {
    my ($self) = @_;

    my @lines = map {"$_"} @{$self->ramachandrans};
    return @lines;
}

sub TO_JSON {
    my ($self) = @_;
    return \@{$self->ramachandrans};
}

__PACKAGE__->meta->make_immutable;
1;
