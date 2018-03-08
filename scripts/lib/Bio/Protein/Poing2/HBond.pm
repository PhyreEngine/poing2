package Bio::Protein::Poing2::HBond;
use strict;
use warnings;
use utf8;
use Carp;
use Bio::Protein::Poing2::Atom;
use namespace::autoclean;
use Moose;

=head1 NAME

Bio::Protein::Poing2::HBond - Hydrogen bond between two atoms

=head1 SYNOPSIS

    #H bond between $donor and $acceptor residues with a distance of 3.3A.  The
    #angle between the donor N atom and the acceptor O and C atoms is given by
    #the "noc" field. The angle between the O atom of the acceptor, the N atom
    #of the donor and the C atom of the residue before the donor is given by
    #the "onc" field.

    #The residue before the donor atom is the C<prev> residue. If this is not
    #supplied, then the C<onc_angle> method will throw an exception.

    my $hbond = Bio::Protein::Poing2::LinearSpring->new(
        donor => $donor_residue,
        acceptor => $acceptor_residue,
        prev => $prev_residue,
        distance => 3.3, #Angstrom
        noc => 107, #degrees
        onc => 125, #degrees
    );
    #Get a Bio::Protein::Poing2::LinearSpring
    $hbond->linear;
    #Get Bio::Protein::Poing2::BondAngle objects
    $hbond->noc_angle;
    $hbond->onc_angle;

=head1 METHODS

=over

=item C<new(%args)>: Instantiate a new C<HBond> object.

Requires the following arguments:

=over

=item C<donor> (L<Bio::Protein::Poing2::Residue): donor residue.

=item C<acceptor> (L<Bio::Protein::Poing2::Residue): acceptor residue.

=item C<prev> (Optional L<Bio::Protein::Poing2::Residue): Residue before the
donor residue.

=item C<distance>: N-O distance in Angstroms.

=item C<noc>: N..O=C angle in degrees.

=item C<onc>: O..N-C angle in degrees.

=back

=cut

has donor    => (is => 'ro', required => 1, isa => 'Bio::Protein::Poing2::Residue');
has acceptor => (is => 'ro', required => 1, isa => 'Bio::Protein::Poing2::Residue');
has prev     => (is => 'ro', required => 0, isa => 'Maybe[Bio::Protein::Poing2::Residue]');

has distance => (is => 'ro', isa => 'Num',  required => 1);
has onc => (is => 'ro', isa => 'Num',  required => 1);
has noc => (is => 'ro', isa => 'Num',  required => 1);

=item C<linear>: Get a L<Bio::Protein::Poing2::LinearSpring> representing the
N-O distance.

=cut

sub linear {
    my ($self) = @_;

    return Bio::Protein::Poing2::LinearSpring->new(
        atom_1 => $self->donor->atom_by_name('N'),
        atom_2 => $self->acceptor->atom_by_name('O'),
        distance => $self->distance,
    );
}

=item C<noc_angle>: Get a L<Bio::Protein::Poing2::BondAngle> for the N..O=C
angle.

=cut

sub noc_angle {
    my ($self) = @_;
    return Bio::Protein::Poing2::BondAngle->new(
        atoms => [
            $self->donor->atom_by_name('N'),
            $self->acceptor->atom_by_name('O'),
            $self->acceptor->atom_by_name('C'),
        ],
        angle => $self->noc,
    );
}

=item C<onc_angle>: Get a L<Bio::Protein::Poing2::BondAngle> for the O..N-C
angle.

If the C<prev> residue was not set, this method will die.

=cut

sub onc_angle {
    my ($self) = @_;
    return Bio::Protein::Poing2::BondAngle->new(
        atoms => [
            $self->acceptor->atom_by_name('O'),
            $self->donor->atom_by_name('N'),
            $self->prev->atom_by_name('C'),
        ],
        angle => $self->onc,
    );
}

=back

=cut

__PACKAGE__->meta->make_immutable;
1;
