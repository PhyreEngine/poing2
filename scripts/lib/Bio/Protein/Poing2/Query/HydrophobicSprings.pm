package Bio::Protein::Poing2::Query::HydrophobicSprings;
use Moose::Role;
use Bio::Protein::Poing2::LinearSpring;

=head1 NAME

Bio::Protein::Poing2::Query::HydrophobicSprings - Add springs between hydrophobic residues.

=head1 SYNOPSIS

    #In your package
    use Moose;
    with 'Bio::Protein::Poing2::Query::HydrophobicSprings';
    #Arrayref of Bio::Protein::Poing2::LinearSpring
    my $hydrophobic_springs = $self->hydrophobic_springs;

=head1 DESCRIPTION

This module attempts to force hydrophobic packing by placing springs between
all pairs of hydrophobic residues. We use the Kyre-Doolittle scale here for
hydrophicity. A spring is placed between each pair of residues for which both
the sum and product of the hydrophicity scores for each residue is greater than
zero.

These springs are very weak, and have a fairly arbitrary equilibrium distance
of 8Å. The cutoff is also set to 8Å so that residues must come fairly close to
be attracted to one another.

=cut


requires 'residues';

#Hydrophobicity scale. Positive means more hydrophobic.
our %kyte_doolittle = (
    ALA =>  1.800,
    ARG => -4.500,
    ASN => -3.500,
    ASP => -3.500,
    CYS =>  2.500,
    GLN => -3.500,
    GLU => -3.500,
    GLY => -0.400,
    HIS => -3.200,
    ILE =>  4.500,
    LEU =>  3.800,
    LYS => -3.900,
    MET =>  1.900,
    PHE =>  2.800,
    PRO => -1.600,
    SER => -0.800,
    THR => -0.700,
    TRP => -0.900,
    TYR => -1.300,
    VAL =>  4.200,
);

sub _build_hydrophobic_springs {
    my ($self) = @_;

    my @springs = ();
    for my $i(sort {$a <=> $b} keys %{$self->residues}){
        for my $j(sort {$a <=> $b} keys %{$self->residues}){
            next unless $j < $i;
            my $r_i = $self->residues->{$i};
            my $r_j = $self->residues->{$j};

            my $h_i = $kyte_doolittle{$r_i->threeletter};
            my $h_j = $kyte_doolittle{$r_j->threeletter};

            my $h_sum  = $h_i + $h_j;
            my $h_prod = $h_i * $h_j;

            if($h_sum > 0 && $h_prod > 0){
                my $spring = Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $r_i->atom_by_name($r_i->threeletter),
                    atom_2 => $r_j->atom_by_name($r_j->threeletter),
                    distance => 8,
                    cutoff => 8,
                    constant => 0.005,
                );
                push @springs, $spring;
            }
        }
    }
    return \@springs;
}

1;
