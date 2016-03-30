package Bio::Protein::Poing2::Filter::Aln::Known;
use strict;
use warnings;
use Bio::Protein::Poing2::Data;
use Bio::Protein::Poing2;
use Moose;


=head1 NAME

Bio::Protein::Poing2::Filter::Aln::Known - Remove unknown (UNK/X) residues from an alignment

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Aln::Known->new();
    my $filtered_residues = $filt->filter($aln);

=cut

sub filter {
    my ($self, $aln) = @_;

    for my $res_key(keys %{$aln}){

        my $from = $aln->{$res_key}->{from};
        my $to   = $aln->{$res_key}->{to};

        if($from && $from->threeletter eq 'UNK'){
            $aln->{$res_key}->{from} = undef;
        }
        if($to && $to->threeletter eq 'UNK'){
            $aln->{$res_key}->{to} = undef;
        }
    }
    return $aln;
}

__PACKAGE__->meta->make_immutable;
1;
