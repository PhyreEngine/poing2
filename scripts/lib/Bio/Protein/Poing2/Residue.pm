package Bio::Protein::Poing2::Residue;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Data qw(%CA_SC_len %BB_BB_len);
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::LinearSpring;

use Bio::Protein::Poing2::Class;
use if $^V gt v5.10.1, parent => 'Bio::Protein::Poing2::Class';
use if $^V le v5.10.1, base   => 'Bio::Protein::Poing2::Class';


#Allow overloading to string
use overload q{""} => 'threeletter';

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

has ss_state => (is => 'ro');

=item C<index()>: The index of the residue in the sequence.

=cut

has index => (is => 'ro', required => 1);

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

=item C<add_sidechain()>: Add a sidechain atom to this residue.

=cut

sub add_sidechain {
    my ($self) = @_;
    return if $self->threeletter eq 'GLY';

    push @{$self->atoms}, Bio::Protein::Poing2::Atom->new(
        name    => $self->threeletter,
        residue => $self,
    );
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
}

=back

=cut
1;
