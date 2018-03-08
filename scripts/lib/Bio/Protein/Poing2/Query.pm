package Bio::Protein::Poing2::Query;
use strict;
use warnings;
use Bio::Protein::Poing2::IO::Fasta;
use Bio::Protein::Poing2::IO::Psipred;
use Bio::Protein::Poing2::IO::PDB;
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
has bb_only  => (is => 'ro', isa => 'Bool', default => 0);
has explicit_ss  => (is => 'ro', isa => 'Bool', default => 0);
has hydrophobic_springs  => (is => 'ro', isa => 'Bool', default => 0);

has positions_file => (is => 'ro', isa => 'Maybe[Str]', default => undef);

has backbone_constraints => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_backbone_constraints');
has backbone_springs => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_backbone_springs');
has backbone => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_backbone');
has residues => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_sequence');
has angles   => (is => 'ro', lazy => 1, builder => '_build_angles');
has ramachandran => (is => 'ro', lazy => 1, builder => '_build_ramachandran');
has fourmers => (is => 'ro', lazy => 1, builder => '_init_bb_fourmers');

with 'Bio::Protein::Poing2::Query::HydrophobicSprings';

sub length {
    my ($self) = @_;
    return scalar(keys %{$self->residues});
}

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
    #Add atoms to the residue
    return $residues;
}

sub _build_backbone {
    my ($self) = @_;

    #Add atoms to each residue
    my $residues = $self->residues;
    my $atom_index = 1;
    for my $res_index(sort {$a <=> $b} keys %{$residues}){
        $atom_index = $residues->{$res_index}->init_fine_bb($atom_index);
        if(!$self->bb_only){
            $atom_index = $residues->{$res_index}->init_coarse_sc($atom_index);
        }
    }
    if($self->positions_file){
        $self->_add_positions($residues);
    }
    return $residues;
}

sub _build_backbone_constraints {
    my ($self) = @_;

    #Build the springs for the backbone
    my $springs = $self->_init_bb_springs(
        \%Bio::Protein::Poing2::Data::fine_bb_links,
        \%Bio::Protein::Poing2::Data::BB_BB_len,
    );
    if(!$self->bb_only){
        push @{$springs}, @{$self->_init_bb_springs(
            \%Bio::Protein::Poing2::Data::CA_SC_links,
            \%Bio::Protein::Poing2::Data::CA_SC_len,
        )};
    };
    return $springs;
}

sub _build_backbone_springs {
    my ($self) = @_;

    my $springs = [];
    if($self->ss && $self->explicit_ss){
        push @{$springs}, @{$self->_init_bb_ss_springs};
    }
    if($self->hydrophobic_springs){
        push @{$springs}, @{$self->_build_hydrophobic_springs};
    }

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

=item C<_init_bb_springs($links, $lengths)>:

Initialise backbone springs according to some link and length specification.
See, for example, C<%Bio::Protein::Poing2::Data::fine_bb_links> and
C<%Bio::Protein::Poing2::Data::BB_BB_len>;

=cut

sub _init_bb_springs {
    my ($self, $links, $lengths) = @_;

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
                if($atom_j){
                    push @springs, Bio::Protein::Poing2::LinearSpring->new(
                        atom_1 => $atom_i,
                        atom_2 => $atom_j,
                        distance => $lengths->{$atom_i_name}->{$link->{atom}},
                    );
                }
            }
        }
    }
    return \@springs;
}

sub _init_bb_ss_springs {
    my ($self) = @_;
    my @springs = ();
    for my $res_idx(sort {$a <=> $b} keys %{$self->backbone}){
        my $r1 = $self->residues->{$res_idx};
        my $r2 = $self->residues->{$res_idx + 1};
        my $r3 = $self->residues->{$res_idx + 2};
        my $r4 = $self->residues->{$res_idx + 3};
        my $r5 = $self->residues->{$res_idx + 4};

        if($r1 && $r2 && $r3 && $r4){
            next if !$r1->ss_state || $r1->ss_state ne 'H';
            next if !$r2->ss_state || $r2->ss_state ne 'H';
            next if !$r3->ss_state || $r3->ss_state ne 'H';
            next if !$r4->ss_state || $r4->ss_state ne 'H';
            next if !$r5->ss_state || $r4->ss_state ne 'H';

            my $ca1 = $r1->atom_by_name('CA');
            my $ca2 = $r2->atom_by_name('CA');
            my $ca3 = $r3->atom_by_name('CA');
            my $ca4 = $r4->atom_by_name('CA');
            my $ca5 = $r5->atom_by_name('CA');
            next if !$ca1 || !$ca2 || !$ca3 || !$ca4 || !$ca5;

            push @springs, Bio::Protein::Poing2::LinearSpring->new(
                atom_1 => $ca1,
                atom_2 => $ca5,
                distance => 6.1,

                handedness => "RIGHT",
                inner_atom => $ca2,
                outer_atom => $ca4,
            );
            push @springs, Bio::Protein::Poing2::LinearSpring->new(
                atom_1 => $ca1,
                atom_2 => $ca4,
                distance => 5.0,

                handedness => "RIGHT",
                inner_atom => $ca2,
                outer_atom => $ca3,
            );
        }
    }
    return \@springs;
}

sub _init_bb_fourmers {
    my ($self) = @_;
    my @fourmers = ();

    for my $res_idx(sort {$a <=> $b} keys %{$self->backbone}){
        my $r1 = $self->residues->{$res_idx - 1};
        my $r2 = $self->residues->{$res_idx};
        my $r3 = $self->residues->{$res_idx + 1};

        if($r1 && $r2){
            my @omega_atoms = ();
            for my $link(@Bio::Protein::Poing2::Data::omega_links){
                my $res = $self->residues->{$res_idx + $link->{increment}} // next;
                my $atom = $res->atom_by_name($link->{atom}) // next;
                push @omega_atoms, $atom;
            }
            if(@omega_atoms == 4){
                push @fourmers, Bio::Protein::Poing2::Fourmer->new(
                    atoms => \@omega_atoms,
                    dihedral => 180,
                );
            }

            if(!$self->bb_only){
                my $prev_c = $r1->atom_by_name('C');
                my $n = $r2->atom_by_name('N');
                my $ca = $r2->atom_by_name('CA');
                my $sc = $r2->atom_by_name($r2->threeletter);

                if($prev_c && $n && $ca && $sc){
                    push @fourmers, Bio::Protein::Poing2::Fourmer->new(
                        atoms => [$prev_c, $n, $ca, $sc],
                        dihedral => 180,
                    );
                }
            }
        }

        if($self->explicit_ss && $r1 && $r2 && $r3){

            next if !$r1->ss_state || $r1->ss_state ne 'H';
            next if !$r2->ss_state || $r2->ss_state ne 'H';
            next if !$r3->ss_state || $r3->ss_state ne 'H';

            my @phi_atoms = ();
            my @psi_atoms = ();
            for my $link(@Bio::Protein::Poing2::Data::phi_links){
                my $res = $self->residues->{$res_idx + $link->{increment}} // next;
                my $atom = $res->atom_by_name($link->{atom}) // next;
                push @phi_atoms, $atom;
            }
            for my $link(@Bio::Protein::Poing2::Data::psi_links){
                my $res = $self->residues->{$res_idx + $link->{increment}} // next;
                my $atom = $res->atom_by_name($link->{atom}) // next;
                push @psi_atoms, $atom;
            }

            next if @phi_atoms != 4 || @psi_atoms != 4;

            push @fourmers, Bio::Protein::Poing2::Fourmer->new(
                atoms => \@phi_atoms,
                dihedral => -60,
            );
            push @fourmers, Bio::Protein::Poing2::Fourmer->new(
                atoms => \@psi_atoms,
                dihedral => -45,
            );
        }
    }
    return \@fourmers;
}

=item C<_build_angles($links)>:

Initialise backbone angles according to some link specification. See, for
example, C<%Bio::Protein::Poing2::Data::fine_bb_links>.

=cut

sub _build_angles {
    my ($self) = @_;
    my @angles = ();
    push @angles, @{$self->_build_bb_angles};
    push @angles, @{$self->_build_sc_angles} if !$self->bb_only;
    return \@angles;
}

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

sub _build_sc_angles {
    my ($self) = @_;

    my @angles = ();
    for my $res_idx(sort {$a <=> $b} keys %{$self->backbone}){
        my $res = $self->backbone->{$res_idx};
        if($Bio::Protein::Poing2::Data::SC_angles{$res->threeletter}){
            my $n_atom = $res->atom_by_name('N');
            my $ca_atom = $res->atom_by_name('CA');
            my $sc_atom = $res->atom_by_name($res->threeletter);

            push @angles, Bio::Protein::Poing2::BondAngle->new(
                atoms => [$n_atom, $ca_atom, $sc_atom],
                angle => $Bio::Protein::Poing2::Data::SC_angles{
                    $res->threeletter
                },
            ) if $n_atom && $ca_atom && $sc_atom;
        }
    }
    return \@angles;
}

#Open the PDB file supplied as the "positions_file" argument and read positions
#for each atom from that file. No error checking is done to make sure that all
#atoms have a position set.
sub _add_positions {
    my ($self, $residues) = @_;

    my $pdb_residues = Bio::Protein::Poing2::IO::PDB::read_pdb(
        $self->positions_file);

    for my $res_num(keys %{$pdb_residues}){
        my $q_res = $residues->{$res_num} // next;

        for my $atom(@{$pdb_residues->{$res_num}->atoms}){
            my $q_atom = $q_res->atom_by_name($atom->name) // next;
            $q_atom->coords($atom->coords);
        }
    }
}


__PACKAGE__->meta->make_immutable;
1;
