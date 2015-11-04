package Bio::Protein::Poing2::Template;
use strict;
use warnings;
use Bio::Protein::Poing2::IO::FastaAlignment;
use Bio::Protein::Poing2::IO::PDB;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Fourmer;
use Bio::Protein::Poing2::Ramachandran;
use Bio::Protein::Poing2::Handedness;
use Bio::Protein::Poing2::Filter::Atom::Backbone;
use List::Util qw(min max);
use Moose;

=head1 NAME

Bio::Protein::Poing2::Template - All information from a template

=head1 SYNOPSIS

    use Bio::Protein::Poing2::Template;
    my $template = Bio::Protein::Poing2::Template->new(
        alignment => 'alignment.fasta',
        model     => 'model.pdb',
    );

=head1 DESCRIPTION

When Poing2 is building a model, it can use information from templates. These
templates are built by homology modelling systems, and so we can glean some
information from the alignment and model.

=head2 Constraints

=head3 Pairwise distances

The strongest constraint in Poing2 is the matrix of CA-CA distances found from
each template. These constrain the structure at all ranges, and encode for all
information except handedness.

=head3 Fourmers

If the dihedral angles of a protein can be found, then the protein can be well
approximated. Poing2 uses the dihedral angles of the templates to help
constrain the model under construction. This ensures that handedness remains
correct.

The omega fourmer can always be taken, but the phi/psi constraints are violated
when a residue changes Ramachandran type; that is, when a residue of the
general type changes to a proline or a glycine, and so on. In these cases, we
do not take the fourmers as constraints.

=back

=cut

has alignment => (is => 'ro', isa => 'Str', required => 1);
has model     => (is => 'ro', isa => 'Str', required => 1);
has fourmers  => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_fourmers');
has pairs     => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_pairwise_springs');
has residues  => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_residues');
has aln       => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_alignment');
has handedness => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_handedness');

sub _read_residues {
    my ($self) = @_;

    #We only want the backbone atoms
    my $bb_filter = Bio::Protein::Poing2::Filter::Atom::Backbone->new();

    my $residues = Bio::Protein::Poing2::IO::PDB::read_pdb($self->model);
    $residues = $bb_filter->filter($residues);
    return $residues;
}

sub _build_handedness {
    my ($self) = @_;

    my $residues = $self->residues;
    my $springs = $self->pairs;

    my @hands = ();
    my @res_idx = sort {$a <=> $b} keys %{$residues};
    for my $spring(@{$springs}){
        my $r1 = $spring->atom_1->residue->index;
        my $r2 = $spring->atom_2->residue->index;

        my $sep = int(($r2->index - $r1->index) / 3);
        my $inner = $self->closest_residue($r1 + $sep)->atom_by_name('CA');
        my $outer = $self->closest_residue($r2 - $sep)->atom_by_name('CA');
        next if(!$inner || !$outer);

        my $hand = Bio::Protein::Poing2::Handedness->new(
            atoms => [$spring->atom_1, $inner, $outer, $spring->atom_2],
        );
        push @hands, $hand;
    }
    return \@hands;
}

sub _build_pairwise_springs {
    my ($self) = @_;

    my $pairs = [];

    my $res = $self->residues;

    my $done = 0;
    for my $r1(sort {$a <=> $b} keys %{$res}){
        my $res1 = $res->{$r1};
        my $a1 = $res1->atom_by_name('CA');

        for my $r2(sort {$a <=> $b} keys %{$res}){
            next if $r1 == $r2;
            my $res2 = $res->{$r2};
            my $a2 = $res2->atom_by_name('CA');

            my $pair = Bio::Protein::Poing2::LinearSpring->new(
                atom_1 => $a1,
                atom_2 => $a2,
                distance => abs($a1->coords - $a2->coords),
            );
            push @{$pairs}, $pair;
        }
    }
    return $pairs;
}

sub _read_alignment {
    my ($self) = @_;
    return Bio::Protein::Poing2::IO::FastaAlignment::read_fasta(
        $self->alignment
    );
}

sub _build_fourmers {
    my ($self) = @_;

    #Build fourmers for all omega angles, and for phi/psi angles for which the
    #Ramachandran type has not changed.

    my @fourmers = ();
    my $done = 0;
    my $nres = keys %{$self->residues};
    for my $residue_idx(sort {$a <=> $b} keys %{$self->residues}){

        my $omega = $self->build_fourmer(
            $residue_idx,
            \@Bio::Protein::Poing2::Data::omega_links);
        push @{$self->{omega}}, $omega if defined $omega;

        #check if the residue type has changed (e.g. general -> GLY)
        my $aln = $self->aln->{$residue_idx};
        next if !$aln->{to};

        my $rama_from = Bio::Protein::Poing2::Ramachandran->new(
            residue => $aln->{from}
        );
        my $rama_to   = Bio::Protein::Poing2::Ramachandran->new(
            residue => $aln->{to}
        );
        next if $rama_from->type ne $rama_to->type;

        #See if we can add a fourmer. We can only do so if all the atoms exist
        my $phi = $self->build_fourmer(
            $residue_idx,
            \@Bio::Protein::Poing2::Data::phi_links
        );
        push @{$self->{phi}}, $phi if defined $phi;
        my $psi = $self->build_fourmer(
            $residue_idx,
            \@Bio::Protein::Poing2::Data::psi_links
        );
        push @{$self->{psi}}, $psi if defined $psi;
    }
    return [@{$self->{phi}}, @{$self->{psi}}, @{$self->{omega}}];
}

sub build_fourmer {
    my ($self, $index, $link) = @_;

    #So, $link looks like:
    #[{increment => 0, atom => 'C'}, ...]

    #Get all atoms
    my @atoms = ();
    for my $atom_spec(@{$link}){
        my $residue = $self->residues->{$index + $atom_spec->{increment}};
        next unless $residue;
        my $atom = $residue->atom_by_name($atom_spec->{atom});
        push @atoms, $atom;
    }
    if(@atoms == 4){
        return Bio::Protein::Poing2::Fourmer->new(atoms => \@atoms);
    }else{
        return undef;
    }
}

#Find the closest residue that actually exists in the template
sub closest_residue {
    my ($self, $resi) = @_;

    my $first_resi = min(keys %{$self->residues});
    my $last_resi  = max(keys %{$self->residues});
    for(my $i=0; $resi - $i >= $first_resi || $resi + $i <= $last_resi; $i++){
        #Prefer lower numbers
        my $r1 = $self->residues->{$resi - $i};
        my $r2 = $self->residues->{$resi + $i};
        return $r1 if $r1;
        return $r2 if $r2;
    }
    return undef;
}

__PACKAGE__->meta->make_immutable;
1;
