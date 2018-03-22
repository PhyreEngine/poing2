package Bio::Protein::Poing2::Residue;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data qw(%CA_SC_len %BB_BB_len @AA1);
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::LinearSpring;
use Bio::Protein::Poing2;
use Moose;

#Atoms for a course representation
our %coarse_sidechain = (
    A => [qw(ALA)], C => [qw(CYS)], D => [qw(ASP)], E => [qw(GLU)],
    F => [qw(PHE)], G => [qw(   )], H => [qw(HIS)], I => [qw(ILE)],
    K => [qw(LYS)], L => [qw(LEU)], M => [qw(MET)], N => [qw(ASN)],
    P => [qw(PRO)], Q => [qw(GLN)], R => [qw(ARG)], S => [qw(SER)],
    T => [qw(THR)], V => [qw(VAL)], W => [qw(TRP)], Y => [qw(TYR)],
    Z => [qw(GLX)], B => [qw(ASX)],
);
our %coarse_bb = map {$_ => ['CA']} @AA1;

=head1 NAME

Bio::Protein::Poing2::Residue - Class representing a residue

=head1 SYNOPSIS

    my $res1 = Bio::Protein::Poing2::Residue->new('ALA');
    my $res2 = Bio::Protein::Poing2::Residue->new('A');

    print $res1->oneletter;   #A
    print $res1->threeletter; #ALA
    print "$res1";            #ALA

=head1 METHODS

=over

=item C<new(type => $type)> Instantiate a new class of type C<$type>.

=cut

has type => (is => 'ro', required => 1);

=item C<atoms()>: Returns an arrayref of L<Bio::Protein::Poing2::Atom> objects.

=cut

has atoms => (is => 'rw', default => sub{[]});

=item C<ss_state()>: Returns a single-letter code for the
secondary structure state.

=cut

has ss_state => (is => 'rw', isa => 'Maybe[Str]', default => undef);

=item C<index()>: The index of the residue in the sequence.

=cut

has index => (is => 'rw', required => 1);

=item C<oneletter()>: Returns the one-letter code for this residue.

=cut

sub oneletter {
    my ($self) = @_;
    my $aa = $self->type;
    $aa = $Bio::Protein::Poing2::Data::three2one{$aa} if length $aa == 3;
    return $aa;
}

=item C<threeletter()>: Returns the three-letter code for this residue.

=cut

sub threeletter {
    my ($self) = @_;
    my $aa = $self->type;
    $aa = $Bio::Protein::Poing2::Data::one2three{$aa} if length $aa == 1;
    return $aa;
}

=item C<atom_by_name($name)>: Get atom by name.

=cut

sub atom_by_name {
    my ($self, $name) = @_;
    if(!$self->{atoms_by_name}){
        $self->{atoms_by_name}{$_->name} = $_ for @{$self->atoms};
    }
    return $self->{atoms_by_name}{$name};
}

=item C<add_sidechain([$index])>: Add a sidechain atom to this residue.

If C<$index> is supplied, the added atoms will have the index starting at
C<$index>.

=cut

sub add_sidechain {
    my ($self, $index) = @_;
    $index //= 1;
    return if $self->threeletter eq 'GLY';

    push @{$self->atoms}, Bio::Protein::Poing2::Atom->new(
        name    => $self->threeletter,
        residue => $self,
        index => $index,
    );
}

=item C<add_atom($atom)>: Add C<$atom> to the residue.

=cut

sub add_atom {
    my ($self, $atom) = @_;
    push @{$self->atoms}, $atom;
}

=item C<init_coarse_bb([$start])>: Add coarse backbone atoms.

If C<$start> is supplied, the index of the added atoms begins at C<$start>.

=cut

sub init_coarse_bb {
    my ($self, $start) = @_;
    my $natoms = $self->init_atoms(['CA'], $start);
    return ($start // 1) + $natoms;
}

=item C<init_fine_bb([$start])>: Add fine backbone atoms.

If C<$start> is supplied, the index of the added atoms begins at C<$start>.

Returns the index of the next atom.

=cut

sub init_fine_bb {
    my ($self, $start) = @_;
    my $natoms = $self->init_atoms(
        \@Bio::Protein::Poing2::Data::backbone_order,
        $start);
    return ($start // 1) + $natoms;
}

=item C<init_coarse_sc([$start])>: Add coarse sidechain atoms.

If C<$start> is supplied, the index of the added atoms begins at C<$start>.

=cut

sub init_coarse_sc {
    my ($self, $start) = @_;
    my $natoms = $self->init_atoms(
        $coarse_sidechain{$self->oneletter},
        $start);
    return ($start // 1) + $natoms;
}

=item C<init_atoms($list, [$start]): Add atoms in C<$list>.

Returns the number of atoms added (which may be less than C<@{$list}>, because
duplicate atoms are ignored).

=cut

sub init_atoms {
    my ($self, $atom_names, $start) = @_;
    $start //= 1;

    #Get current atoms so we don't duplicate anything
    my %current = map {$_->name => $_} @{$self->atoms};

    #Add atoms
    my $i = 0;
    for my $name(@{$atom_names}){
        next if $current{$name};

        $self->add_atom(Bio::Protein::Poing2::Atom->new(
            residue => $self,
            name => $name,
            index => $start++,
        ));
        $i++;
    }
    return $i;
}

=item C<internal_springs()>: Get internal springs for this residue.

This includes springs between CA and sidechain atoms, and springs connecting
backbone atoms.

=cut

sub internal_springs {
    my ($self) = @_;

    my @springs = ();
    for my $a1(@{$self->atoms}){
        for my $a2(@{$self->atoms}){
            last if $a1 == $a2;

            if($a1->name eq 'CA' && $CA_SC_len{$a2->name}){
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $a1,
                    atom_2 => $a2,
                    distance => $CA_SC_len{$a2->name}
                );
            }elsif($a2->name eq 'CA' && $CA_SC_len{$a1->name}){
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $a1,
                    atom_2 => $a2,
                    distance => $CA_SC_len{$a1->name}
                );
            }elsif($BB_BB_len{$a1->name} && $BB_BB_len{$a1->name}->{$a2->name}){
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $a1,
                    atom_2 => $a2,
                    distance => $BB_BB_len{$a1->name}->{$a2->name}
                );
            }elsif($BB_BB_len{$a2->name} && $BB_BB_len{$a2->name}->{$a1->name}){
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $a1,
                    atom_2 => $a2,
                    distance => $BB_BB_len{$a2->name}->{$a1->name}
                );
            }
        }
    }
    return \@springs;
}

=item C<string_repr()>: Get a string representation for the config file.

=cut

sub string_repr {
    my ($self) = @_;
    return join q{}, map {$_->string_repr()} @{$self->atoms};
}

use overload '""' => \&threeletter;

__PACKAGE__->meta->make_immutable;
1;
