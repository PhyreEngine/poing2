package Bio::Protein::Poing2;
use strict;
use warnings;
use Bio::Protein::Poing2::LinearSpring;
use Bio::Protein::Poing2::Fourmer;
use Bio::Protein::Poing2::IO::Fasta;
use Bio::Protein::Poing2::IO::PDB;
use Bio::Protein::Poing2::Data qw(%BB_BB_len);
use Moose;

=head1 NAME

Bio::Protein::Poing2 - Generate input files for poing2

=head1 SYNOPSIS

=head1 ATTRIBUTES

=over

=item C<pdbs>

Arrayref of PDB files from which to take constraints.

=cut


has pdbs  => (is => 'ro', default => sub { [] });

=item C<query>

Location of a FASTA file from which to parse the query.

=cut

has query => (is => 'ro', required => 1);

=item C<verbose>

If true, write informative messages to STDERR.

=cut

has verbose => (is => 'rw', default => 0);

=item C<pair_filters>

Arrayref of pair filters.

=cut

has pair_filters => (is => 'ro', default => sub { [] });

=item C<atom_filters>

Arrayref of atom filters.

=cut

has atom_filters => (is => 'ro', default => sub { [] });

=back

=head1 METHODS

=over

=item C<residues()>

Get a hashref of array residues mapping residue index to a Residue object.

=cut

sub residues {
    my ($self) = @_;
    $self->{residues} ||= Bio::Protein::Poing2::IO::Fasta::read_fasta(
        $self->query
    );
    return $self->{residues};
}

=item C<template_residues()>

Get a hashref of residues. The keys of the hashref are the template names, and
the values are a hashref of residues:

    {
        template1 => {
            1 => Residue object,
            2 => Residue object,
            ...
        },
        template2 => {...}
        ...
    }

=cut

sub template_residues {
    my ($self) = @_;
    return $self->{template_residues} if $self->{template_residues};

    my %template_res = ();
    for my $template(@{$self->pdbs}){
        my $res = Bio::Protein::Poing2::IO::PDB::read_pdb($template);

        for my $filter(@{$self->atom_filters}){
            $res = $filter->filter($res);
        }
        $template_res{$template} = $res;
    }
    $self->{template_residues} = \%template_res;
    return $self->{template_residues};
}

=item C<set_positions()>

Set positions of atoms from templates.

=cut

sub set_positions {
    my ($self) = @_;

    my $residues = $self->residues;
    for my $t1(@{$self->pdbs}){
        my $res_list = $self->template_residues->{$t1};

        for my $r1(sort {$a <=> $b} keys %{$res_list}){
            my $template_res = $res_list->{$r1};
            my $set_res = $residues->{$r1};
            next unless defined $set_res;

            for my $template_atom(@{$template_res->atoms}){
                #Atom to set coords of
                my $to_set = $set_res->atom_by_name($template_atom->name);
                $to_set->coords($template_atom->coords) if $to_set;
            }
        }
    }
}

=item C<pairs()>

Get an arrayref of residue pairs, represented by
L<Bio::Protein::Poing2::LinearSpring> objects. These will be generated between
each residue pair in the templates, filtered by the contents of
C<pair_filters>.

=cut

sub pairs {
    my ($self) = @_;
    return $self->{pairs} if $self->{pairs};

    my $pairs = [];
    for my $t1(@{$self->pdbs}){

        my $res = $self->template_residues->{$t1};
        my $done = 0;
        for my $r1(sort {$a <=> $b} keys %{$res}){
            print STDERR "\rBuilding pairs for $t1: ",
                sprintf("%.1f", ((++$done) / keys %{$res}) * 100),
                '%'
                if $self->verbose;

            for my $r2(sort {$a <=> $b} keys %{$res}){
                next if $r1 == $r2;

                for my $a1(@{$res->{$r1}->atoms}){
                    for my $a2(@{$res->{$r2}->atoms}){
                        push @{$pairs}, Bio::Protein::Poing2::LinearSpring->new(
                            atom_1 => $a1,
                            atom_2 => $a2,
                            distance => abs($a1->coords - $a2->coords),
                        );
                    }
                }
            }
        }
        print STDERR "\n" if $self->verbose;
    }

    for my $filter(@{$self->pair_filters}){
        $pairs = $filter->filter($pairs);
    }
    $self->{pairs} = $pairs;
    return $self->{pairs};
}

=item C<fourmers()>

Get an arrayref of fourmers.

=cut

sub fourmers {
    my ($self) = @_;
    return $self->{fourmers} if defined $self->{fourmers};

    my @fourmers = ();
    push @fourmers, $self->all_fourmers(\@Bio::Protein::Poing2::Data::phi_links);
    push @fourmers, $self->all_fourmers(\@Bio::Protein::Poing2::Data::psi_links);
    $self->{fourmers} = \@fourmers;

    return $self->{fourmers};
}

sub all_fourmers {
    my ($self, $spec) = @_;

    my @fourmers = ();

    #Build fourmer for Phi/Psi angles
    for my $template_id(keys %{$self->template_residues}){
        my $template = $self->template_residues->{$template_id};


        my $done = 0;
        my $nres = keys %{$template};
        for my $residue_idx(keys %{$template}){
            if($self->verbose){
                printf STDERR "\rBuilding fourmers for %s: %2.1f%%",
                    $template_id,
                    (++$done) / $nres * 100;
            }
            #See if we can add a phi fourmer.
            my $fourmer = $self->build_fourmer($template, $residue_idx, $spec);
            push @fourmers, $fourmer if defined $fourmer;
        }
        print STDERR "\n" if $self->verbose;
    }
    return @fourmers;
}

#Look at each residue and those residues behind it to see if we can add
#a fourmer for this. The fourmer is specified by four hashrefs, giving
#the atom and the number of residues behind the current residue. For example,
#
#   {increment => -1, atom => 'N'}
#
#means atom N of the previous residue. So to build all fourmers, we
#scan through all residues in each template, then scan through the
#specs for the fourmers we want to add. If all the specified atoms
#exist, then we add the fourmer to the list of possible fourmers.

sub build_fourmer {
    my ($self, $template, $residue_idx, $spec_list) = @_;

    my @atoms = ();
    for my $spec(@{$spec_list}){
        #Get the specified residue and atoms
        my $other_res = $template->{$residue_idx + $spec->{increment}};
        last if not defined $other_res;
        my %other_atoms = map {$_->name => $_} @{$other_res->atoms};

        my $spec_atom = $other_atoms{$spec->{atom}};
        push @atoms, $spec_atom if defined $spec_atom;
    }
    if(@atoms == 4){
        my $fourmer = Bio::Protein::Poing2::Fourmer->new(
            atoms => \@atoms
        );
        return $fourmer;
    }
    return undef;
}

sub string_repr {
    my ($self) = @_;

    my @lines;
    push @lines, "[PDB]\n";
    for my $i(sort {$a <=> $b} keys %{$self->residues}){
        push @lines, $self->residues->{$i}->string_repr;
    }

    push @lines, "[Linear]\n";
    for my $pair(@{$self->pairs}){
        push @lines, $pair->string_repr;
    }

    push @lines, "[Torsion]\n";
    for my $torsion(@{$self->fourmers}){
        push @lines, $torsion->string_repr;
    }
    return join q{}, @lines;
}

=item C<renumber_atoms()>:

Fix atom numbers to be consecutive.

=cut

sub renumber_atoms {
    my ($self) = @_;

    my $i = 1;
    for(sort {$a <=> $b} keys %{$self->residues}){
        my $res = $self->residues->{$_};
        $_->index($i++) for @{$res->atoms};
    }
}

=item C<init_coarse_bb()>:

Initialise residues with coarse backbone representation.

=cut

sub init_coarse_bb {
    my ($self) = @_;

    #Add atoms to BB
    $_->init_coarse_bb for values %{$self->residues};

    #Add springs
    $self->init_bb_springs(\%Bio::Protein::Poing2::Data::coarse_bb_links);
}

=item C<init_fine_bb()>:

Initialise residues with fine backbone representation.

=cut

sub init_fine_bb {
    my ($self) = @_;

    #Add atoms to BB
    $_->init_fine_bb for values %{$self->residues};

    #Add springs
    $self->init_bb_springs(\%Bio::Protein::Poing2::Data::fine_bb_links);
}

=item C<init_bb_springs($links)>:

Initialise backbone springs according to some link specification. See, for
example, C<%Bio::Protein::Poing2::Data::fine_bb_links>.

=cut

sub init_bb_springs {
    my ($self, $links) = @_;

    for my $res_i_idx(sort {$a <=> $b} keys %{$self->residues}){
        my $res_i = $self->residues->{$res_i_idx};

        # $atom_i will be the atom name in $res_i
        for my $atom_i_name(keys %{$links}){
            for my $link(@{$links->{$atom_i_name}}){
                my $atom_i = $res_i->atom_by_name($atom_i_name);

                #Get connected residue
                my $res_j = $self->residues->{$res_i_idx + $link->{increment}};
                #Residue might not exist
                next if !$res_j;

                #Add spring between atoms i and j
                my $atom_j = $res_j->atom_by_name($link->{atom});
                push @{$self->pairs}, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $atom_i,
                    atom_2 => $atom_j,
                    distance => $Bio::Protein::Poing2::BB_BB_len{$atom_i_name}->{$link->{atom}},
                );
            }
        }
    }
}

=item C<init_coarse_sc()>:

Initialise residues with coarse sidechain representation.

=cut

sub init_coarse_sc {
    my ($self) = @_;
    for(values %{$self->residues}){
        $_->init_coarse_sc;

        my $ca_atom = $_->atom_by_name('CA');
        my $sc_atom = $_->atom_by_name($_->threeletter);
        if($sc_atom && $ca_atom){
            push @{$self->pairs}, Bio::Protein::Poing2::LinearSpring->new(
                atom_1 => $ca_atom,
                atom_2 => $sc_atom,
                distance => $Bio::Protein::Poing2::Data::CA_SC_len{$_->threeletter},
            );
        }

    }
}

__PACKAGE__->meta->make_immutable;
1;