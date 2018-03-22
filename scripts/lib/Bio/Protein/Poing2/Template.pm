package Bio::Protein::Poing2::Template;
use strict;
use warnings;
use Bio::Protein::Poing2::IO::FastaAlignment;
use Bio::Protein::Poing2::IO::PDB;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::HBond;
use Bio::Protein::Poing2::Fourmer;
use Bio::Protein::Poing2::Ramachandran;
use Bio::Protein::Poing2::Filter::Atom::Backbone;
use Bio::Protein::Poing2::Filter::Residue::Known;
use Bio::Protein::Poing2::Filter::Aln::Known;
use List::Util qw(min max);
use List::BinarySearch qw(binsearch_pos);
use Carp;
use autodie;
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
information except angular information.

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

=head1 ATTRIBUTES

The constructor of this class takes the following attributes:

=over

=item C<alignment> (Mandatory string)

Path of a FASTA alignment file containing the alignment of the query to this
template.

=item C<model> (Mandatory string)

Path of a PDB file containing the homology model built for this template.

=item C<query> (Optional L<Bio::Protein::Poing2::Query> object)

If supplied, all the atom IDs of this template will be the same as those in the
query. If atoms are present in the template that are not present in the query,
they will be discarded.

=item C<add_hbonds> (Boolean. Default: false)

Should we run L<stride|http://webclu.bio.wzw.tum.de/stride/> to find hydrogen
bonds and add springs representing these bonds? If the C<stride> executable is
not in the system path, an exception will be thrown.

=back

=cut

has alignment => (is => 'ro', isa => 'Str', required => 1);
has model     => (is => 'ro', isa => 'Str', required => 1);
has add_hbonds=> (is => 'ro', isa => 'Bool', default => 0);
has query     => (is => 'ro', isa => 'Maybe[Bio::Protein::Poing2::Query]', required => 0);
has fourmers  => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_fourmers');
has pairs     => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_pairwise_springs');
has residues  => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_residues');
has hbonds    => (is => 'ro', lazy => 1, init_arg => undef, builder => '_build_hbonds');
has aln       => (is => 'ro', lazy => 1, init_arg => undef, builder => '_read_alignment');

sub _read_residues {
    my ($self) = @_;

    #We only want the backbone atoms
    my $bb_filter = Bio::Protein::Poing2::Filter::Atom::Backbone->new();

    #Discard any unknown residues
    my $known_filter = Bio::Protein::Poing2::Filter::Residue::Known->new();

    my $residues = Bio::Protein::Poing2::IO::PDB::read_pdb($self->model);
    $residues = $bb_filter->filter($residues);
    $residues = $known_filter->filter($residues);

    #If we have a query, remap the residues.
    $residues = $self->_remap_atoms($residues) if $self->query;

    return $residues;
}

sub sorted_residues {
    my ($self) = @_;
    $self->{sorted_residues} ||= [sort {$a <=> $b} keys %{$self->residues}];
    return $self->{sorted_residues};
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
            next if $r1 >= $r2;
            my $res2 = $res->{$r2};
            my $a2 = $res2->atom_by_name('CA');

            #Inner and outer atoms for handedness determination
            my ($inner, $outer);

            if(abs($r2 - $r1) > 4){
                my $sep = int(($r2 - $r1) / 4);
                $inner = $self->closest_residue(
                    $r1 + $sep,
                    round => 'up',
                )->atom_by_name('CA');
                $outer = $self->closest_residue(
                    $r2 - $sep,
                    round => 'down',
                )->atom_by_name('CA');

                if($inner == $a1 || $outer == $a2 || $inner == $outer){
                    $inner = undef;
                    $outer = undef;
                }
            }

            my $pair = Bio::Protein::Poing2::LinearSpring->new(
                atom_1 => $a1,
                atom_2 => $a2,
                inner_atom => $inner,
                outer_atom => $outer,
                distance => abs($a1->coords - $a2->coords),
            );
            push @{$pairs}, $pair;
        }
    }

    #Add H bonds if we were told to
    if($self->add_hbonds){
        my $hbonds = $self->hbonds;
        for my $hbond(@{$hbonds}){
            push @{$pairs}, $hbond->linear;
        }
    }
    return $pairs;
}

sub _build_hbonds {
    my ($self) = @_;

    my @hbonds = ();

    #Run stride:
    my @cmd = ('stride', '-h', $self->model);
    open my $stride_in, q{-|}, @cmd;
    while(my $ln = <$stride_in>){
        next unless $ln =~ /^(?:ACC|DNR)/;
        #Stride H bond lines look like this:
        #REM  |--Residue 1--|     |--Residue 2--|  N-O N..O=C O..N-C     A1     A2  ~~~~
        #ACC  SER A    3    1 ->  TRP A    7    5  3.0  156.0  115.2   37.3   51.8  ~~~~
        #DNR  THR A    4    2 ->  ALA A   30   28  3.0  154.6  118.6   56.8   51.0  ~~~~

        #The first record says that residue 3 is the acceptor (i.e. the O side)
        #and residue 7 is the donor (i.e. the NH side). The second record has
        #the residues in reverse order.

        #The N-O column gives the distance between the N and O atoms (which we
        #will use to add a spring). The "N..O=C" column gives the angle between
        #the N atom of the donor, and the O and C atoms of the acceptor. The
        #"O..N-C" column gives the angle between the O atom of the acceptor,
        #the N atom of the donor and the C atom of the residue before the donor.

        #We will add a spring between the N and O atoms, and add angle
        #constraints for the NOC and ONC atoms.

        my $type = substr $ln, 0, 3;
        my $index1 = substr $ln, 10, 5;
        my $index2 = substr $ln, 30, 5;
        my $no = substr $ln, 40, 5;
        my $noc = substr $ln, 45, 7;
        my $onc = substr $ln, 52, 7;

        $index1 =~ s/ //g;
        $index2 =~ s/ //g;
        $no =~ s/ //g;
        $noc =~ s/ //g;
        $onc =~ s/ //g;

        my $donor_idx = $type eq 'DNR' ? $index1 : $index2;
        my $acc_idx   = $type eq 'DNR' ? $index2 : $index1;

        my $donor = $self->residues->{$donor_idx};
        my $acceptor = $self->residues->{$acc_idx};
        my $prev = $self->residues->{$donor_idx - 1};
        next unless $donor && $acceptor;

        my $n_atm = $donor->atom_by_name('N');
        my $o_atm = $acceptor->atom_by_name('O');
        next unless $n_atm && $o_atm;

        my $hbond = Bio::Protein::Poing2::HBond->new(
            donor => $donor,
            acceptor => $acceptor,
            distance => $no,
            noc => $noc,
            onc => $onc,
            prev => $prev,
        );
        push @hbonds, $hbond;
    }
    close $stride_in;

    return \@hbonds;
}

sub _read_alignment {
    my ($self) = @_;
    my $filter = Bio::Protein::Poing2::Filter::Aln::Known->new();
    my $aln = Bio::Protein::Poing2::IO::FastaAlignment::read_fasta(
        $self->alignment
    );
    $aln = $filter->filter($aln);
    return $aln;
}

sub _build_fourmers {
    my ($self) = @_;

    #Build fourmers for all omega angles, and for phi/psi angles for which the
    #Ramachandran type has not changed.

    $self->{phi}   ||= [];
    $self->{psi}   ||= [];
    $self->{omega} ||= [];

    my @fourmers = ();
    my $done = 0;
    my $nres = keys %{$self->residues};
    # Offset between the model and alignment residue indices. This is
    # increased when a gap in the query is aligned to a residue in the
    # template.
    my $query_offset = 0;
    for my $residue_idx(sort {$a <=> $b} keys %{$self->residues}){

        my $omega = $self->build_fourmer(
            $residue_idx,
            \@Bio::Protein::Poing2::Data::omega_links);
        push @{$self->{omega}}, $omega if defined $omega;

        #check if the residue type has changed (e.g. general -> GLY)
        my $aln = $self->aln->{$residue_idx - $query_offset};

        # If the template residue is aligned to a gap in the query, continue
        if(not defined $aln->{from}) {
            $query_offset++;
            next;
        }

        # Raise an error if the residue type of the model is different to the
        # residue type of this alignment pair.
        if($self->residues->{$residue_idx}->oneletter
               ne  $aln->{from}->oneletter) {
           my $model_type = $self->residues->{$residue_idx}->threeletter;
           my $aln_type = $aln->{from}->threeletter;
           die "Model residue $residue_idx is $model_type, but the "
                . "alignment shows the query sequence as having"
                . " a $aln_type at that position.\n";
        }

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
        #Atom may not exist even if the residue does (CA traces, etc)
        next unless $atom;
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
    my ($self, $resi, %args) = @_;

    my $round = $args{round} || 'up';
    if($self->{closest}{$round}{$resi}){
        return $self->{closest}{$round}{$resi};
    }

    my $residues = $self->sorted_residues;

    my $pos = binsearch_pos {$a <=> $b} $resi, @{$residues};
    if($args{round}){
        if($args{round} eq 'down'){
            $pos = $pos - 1 if $residues->[$pos - 1];
        }
    }
    my $res = $self->residues->{$residues->[$pos]};
    $self->{closest}{$round}{$resi} = $res;
    return $res;
}

#Remap atom IDs to match those of a query. Expects a ::Query object

sub _remap_atoms {
    my ($self, $residues) = @_;

    for my $resi(keys %{$residues}){
        my @replacement_atoms = ();

        for my $atom(@{$residues->{$resi}->atoms}){
            #Check if the atom exists in the query
            my $qres = $self->query->residues->{$resi};
            if($qres){
                my $qatom = $qres->atom_by_name($atom->name);
                $atom->index($qatom->index);
                push @replacement_atoms, $atom;
            }
        }
        $residues->{$resi}->atoms(\@replacement_atoms);
    }
    return $residues;
}

__PACKAGE__->meta->make_immutable;
1;
