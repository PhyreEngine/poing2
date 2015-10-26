package Bio::Protein::Poing2::Query;
use strict;
use warnings;
use Bio::Protein::Poing2::IO::Fasta;
use Bio::Protein::Poing2::IO::Psipred;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Fourmer;
use Bio::Protein::Poing2::BondAngle;
use Bio::Protein::Poing2::Ramachandran::List;
use Moose;


=head1 NAME

Bio::Protein::Poing2::Query - Poing2 representation of a query protein

=head1 SYNOPSIS

    use Bio::Protein::Poing2::Query;
    my $query = Bio::Protein::Poing2::Query->new(
        sequence => 'query.fasta',
        ss => 'psipred.ss2',
    );

=head1 DESCRIPTION

This class contains a representation of the query protein. For the query
protein, we know (or can estimate) backbone constraints, bond angles and
Ramachandran type. We can also make use of secondary structure prediction.

=cut

has sequence => (is => 'ro', isa => 'Str', required => 1);
has ss       => (is => 'ro', isa => 'Maybe[Str]');

has backbone_springs => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_backbone_springs');
has backbone => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_backbone');
has residues => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_sequence');
has angles   => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_bb_angles');
has ramachandran => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_ramachandran');

sub _read_sequence {
    my ($self) = @_;

    #Read sequence from fasta file
    my $residues = Bio::Protein::Poing2::IO::Fasta::read_fasta(
        $self->sequence
    );
    #Now read ss states, if supplied
    if($self->ss){
        my $ss_res = Bio::Protein::Poing2::IO::Psipred::read_psipred(
            $self->ss
        );
        $residues->{$_}->ss_state($ss_res->{$_}->ss_state) for keys %{$ss_res};
    }
    return $residues;
}

sub _build_backbone {
    my ($self) = @_;

    #Add atoms to each residue
    my $residues = $self->residues;
    my $atom_index = 1;
    for my $res_index(sort {$a <=> $b} keys %{$residues}){
        $atom_index = $residues->{$res_index}->init_fine_bb($atom_index);
    }
    return $residues;
}

sub _build_backbone_springs {
    my ($self) = @_;

    #Build the springs for the backbone
    my $springs = $self->_init_bb_springs(
        \%Bio::Protein::Poing2::Data::fine_bb_links
    );
    return $springs;
}

=item C<_build_ramachandran()>:

Calculate Ramachandran regions for this query.

=cut

sub _build_ramachandran {
    my ($self) = @_;
    return Bio::Protein::Poing2::Ramachandran::List->new(
        residues => $self->residues,
    );
}

=item C<_init_bb_springs($links)>:

Initialise backbone springs according to some link specification. See, for
example, C<%Bio::Protein::Poing2::Data::fine_bb_links>.

=cut

sub _init_bb_springs {
    my ($self, $links) = @_;

    my @springs = ();
    for my $res_i_idx(sort {$a <=> $b} keys %{$self->backbone}){
        my $res_i = $self->residues->{$res_i_idx};

        # $atom_i will be the atom name in $res_i
        for my $atom_i_name(keys %{$links}){
            for my $link(@{$links->{$atom_i_name}}){
                my $atom_i = $res_i->atom_by_name($atom_i_name);

                #Get connected residue
                my $res_j = $self->backbone->{$res_i_idx + $link->{increment}};
                #Residue might not exist
                next if !$res_j;

                #Add spring between atoms i and j
                my $atom_j = $res_j->atom_by_name($link->{atom});
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $atom_i,
                    atom_2 => $atom_j,
                    distance => $Bio::Protein::Poing2::Data::BB_BB_len{$atom_i_name}->{$link->{atom}},
                );
            }
        }
    }
    return \@springs;
}

=item C<init_bb_angles($links)>:

Initialise backbone angles according to some link specification. See, for
example, C<%Bio::Protein::Poing2::Data::fine_bb_links>.

=cut

sub _build_bb_angles {
    my ($self) = @_;

    my $links = \%Bio::Protein::Poing2::Data::fine_bb_angles;

    my @angles = ();
    for my $res_i_idx(sort {$a <=> $b} keys %{$self->backbone}){
        my $res_i = $self->backbone->{$res_i_idx};

        #$atom_i_name will be the atom name in $res_i
        for my $atom_i_name(keys %{$links}){
            my $atom_i = $res_i->atom_by_name($atom_i_name);

            for my $link(@{$links->{$atom_i_name}}){

                my $atom_j_name = $link->{atoms}[0]{atom};
                my $atom_k_name = $link->{atoms}[1]{atom};

                my $atom_j_increment = $link->{atoms}[0]{increment};
                my $atom_k_increment = $link->{atoms}[1]{increment};

                my $angle = $link->{angle};

                my $res_j = $self->backbone->{$res_i_idx + $atom_j_increment};
                my $res_k = $self->backbone->{$res_i_idx + $atom_k_increment};

                #Residue might not exist
                next if !$res_j || !$res_k;

                #Add spring between atoms i and j
                my $atom_j = $res_j->atom_by_name($atom_j_name);
                my $atom_k = $res_k->atom_by_name($atom_k_name);

                push @angles, Bio::Protein::Poing2::BondAngle->new(
                    atoms => [$atom_i, $atom_j, $atom_k],
                    angle => $link->{angle},
                );
            }
        }
    }
    return \@angles;
}

__PACKAGE__->meta->make_immutable;
1;
