package Bio::Protein::Poing2::Data;
use strict;
use warnings;
use utf8;

#Allow import of symbols
require Exporter;
use base 'Exporter';

=head1 NAME

Bio::Protein::Poing2::Data - Various constants for Poing2

=head1 SYNOPSIS

    use Bio::Protein::Poing2::Data;
    print $Bio::Protein::Poing2::Data::bond_lens{ALA};

=head1 DESCRIPTION

This module amalgamates a lot of useful constants for Poing2, including various
bond length parameters.

=head1 PACKAGE VARIABLES

=over

=item C<$PI>: π, to machine precision.

=cut

our $PI = 4 * atan2(1, 1);

=item C<@AA1>: Array of single-letter AA codes.

=cut

my @AA1 = qw(A C D E F G H I K L M N P Q R S T V W Y Z B);


=item C<@AA3>: Array of three-letter AA codes;

=cut

my @AA3 = qw(ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR
             VAL TRP TYR GLX ASX);

=item C<%one2three>: Hash of single-letter AA codes to three-letter codes.

=cut

our %one2three = ();
@one2three{@AA1} = @AA3;

=item C<%three2one>: Hash of three-letter AA codes to one-letter codes.

=cut

our %three2one = ();
@three2one{@AA3} = @AA1;

=item C<%CA_SC_len>: CA to sidechain centre of mass mappings, indexed by
three-letter AA code.

=cut

our %CA_SC_len = (
    ALA => 1.52370477561,
    CYS => 2.06811484067,
    ASP => 2.47099367559,
    GLU => 3.09553537293,
    PHE => 3.40440459167,
    HIS => 3.14743193894,
    ILE => 2.31086434835,
    LYS => 3.4901931337 ,
    LEU => 2.60580156453,
    MET => 2.95965999222,
    ASN => 2.47427117216,
    PRO => 1.86904639926,
    GLN => 3.08907392747,
    ARG => 4.11059670604,
    SER => 1.89881631226,
    THR => 1.93550699219,
    VAL => 1.95268816215,
    TRP => 3.86897045069,
    TYR => 3.40900993765,
);

=item B<%backbone>: Atoms in the backbone.

=cut

our %backbone = (
    C  => 2,
    CA => 1,
    N  => 1,
    O  => 1,
);

=item C<%BB_BB_len>: Backbone-backbone distances, stored as a hashref of
distances.

For example, the C-N distance is given by C<$BB_BB_len{C}->{N}> and the C-O
distance is given by C<$BB_BB_len{C}->{O}>. The atoms indexed here are C, N, O
and CA.

=cut

our %BB_BB_len = (
    C  => {N  => 1.3298236446806,  O  => 1.23203299355537},
    N  => {CA => 1.45999906298163                         },
    CA => {C  => 1.53135743668343, CA => 3.8              },
);

=item C<@phi_links>: Links used when calculating the φ dihedral angle (TODO:
reference build_dihedral_sets).

=cut

our @phi_links = (
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
   {increment =>  0, atom => 'CA'},
   {increment =>  0, atom => 'C'},
);

=item C<@psi_links>: Links used when calculating the ψ dihedral angle (TODO:
reference build_dihedral_sets).

=cut

our @psi_links = (
   {increment => -1, atom => 'N'},
   {increment => -1, atom => 'CA'},
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
);

=item C<@omega_links>: Links used when calculating the ω dihedral angle(TODO:
reference build_dihedral_sets).

=cut

our @omega_links = (
   {increment => -1, atom => 'CA'},
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
   {increment =>  0, atom => 'CA'},
);

=item C<@CA_links>: Links used when calculating the torsion angle between the
CA atoms of four residues. (TODO: reference build_dihedral_sets).

=cut

our @CA_links = (
    {increment => -3, atom => 'CA'},
    {increment => -2, atom => 'CA'},
    {increment => -1, atom => 'CA'},
    {increment =>  0, atom => 'CA'},
);

=back

=cut
1;
