package Bio::Protein::Poing2::Filter::Residue::Known;
use strict;
use warnings;
use Bio::Protein::Poing2::Data;
use Bio::Protein::Poing2;
use Moose;


=head1 NAME

Bio::Protein::Poing2::Filter::Residue::Known - Remove unknown (UNK/X) residues

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Residue::Known->new();
    my $filtered_residues = $filt->filter($hashref_of_residues);

=cut

sub filter {
    my ($self, $residues) = @_;
    my %new_res = ();
    for my $res_key(keys %{$residues}){
        if($residues->{$res_key}->threeletter ne 'UNK'){
            $new_res{$res_key} = $residues->{$res_key};
        }
    }
    return \%new_res
}

__PACKAGE__->meta->make_immutable;
1;
