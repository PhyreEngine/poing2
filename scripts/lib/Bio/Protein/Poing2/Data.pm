package Bio::Protein::Poing2::Data;
use strict;
use warnings;
use utf8;

#Allow import of symbols
require Exporter;
use base 'Exporter';
our @EXPORT_OK = qw(%CA_SC_len %BB_BB_len $PI @AA1 @AA3 %fine_bb_links);

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

my @AA1 = qw(A C D E F G H I K L M N P Q R S T V W Y Z B X);


=item C<@AA3>: Array of three-letter AA codes;

=cut

my @AA3 = qw(ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR
             VAL TRP TYR GLX ASX UNK);

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
    CA => {
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
    },
);

=item C<%CA_SC_links>: Links between CA atom and sidechain.

=cut

our %CA_SC_links = (
    CA => [map { {increment => 0, atom => $_} } keys %{$CA_SC_len{CA}}],
);


=item B<%backbone>: Atoms in the backbone.

=cut

our %backbone = (
    C  => 2,
    CA => 1,
    N  => 1,
    O  => 1,
);

=item B<%SC_angles>: Sidechain angles. Only proline is really special.

=cut

our %SC_angles = (
    ALA => 110.2205,
    ARG => 110.5078,
    ASN => 110.3937,
    ASP => 110.3816,
    CYS => 110.4102,
    GLN => 110.4914,
    GLU => 110.4259,
    GLY => 110.1317,
    HIS => 110.4651,
    ILE => 111.1666,
    LEU => 110.1939,
    LYS => 110.4617,
    MET => 110.5773,
    PHE => 110.5342,
    PRO => 103.6294,
    SER => 110.3574,
    THR => 111.1417,
    TRP => 110.3800,
    TYR => 110.4236,
    VAL => 111.2750,
);


=item B<@backbone_order>: Backbone atoms in order.

=cut

our @backbone_order = qw(N CA C O);

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

=item C<%fine_bb_links>: Links between atoms in a fine backbone.

For example, to initialise springs between atoms of a set of residues:

    for my $res_i_idx(keys %residues){
        my $res_i = $residues{$res_i_idx};

        # $atom_i will be the atom name in $res_i
        for my $atom_i_name(keys %Bio::Protein::Poing2::fine_bb_links){
            for my $link(@{$Bio::Protein::Poing2::fine_bb_links{$atom_i_name}}){
                my $atom_i = $res_i->atom_by_name($atom_i_name);

                #Get connected residue
                my $res_j = $residues{$res_i_idx + $link->{increment}};

                #Residue might not exist
                next if !$res_j;

                #Add spring between atoms i and j
                my $atom_j = $res_j->atom_by_name($link->{atom});
                push @springs, Bio::Protein::Poing2::LinearSpring->new(
                    atom_1 => $atom_i,
                    atom_2 => $atom_2,
                    distance => $Bio::Protein::Poing2::BB_BB_len{$atom_i_name}->{$link->{atom}};
                );
            }
        }
    }

=cut

our %fine_bb_links = (
    C  => [{increment => +1, atom => 'N'}, {increment => 0, atom => 'O'}],
    N  => [{increment =>  0, atom => 'CA'}],
    CA => [{increment =>  0, atom => 'C'}],
);

=item C<%coarse_bb_links>: Coarse backbone links. See C<$fine_bb_links>.

=cut

our %coarse_bb_links = (
    CA => [{increment => +1, atom => 'CA'}],
);

our %fine_bb_angles = (
    C => [
        {
            atoms => [
                {increment => +1, atom => 'N'},
                {increment => +1, atom => 'CA'}
            ],
            angle => 121.7,
        },
    ],
    CA => [
        {
            atoms => [
                {increment =>  0, atom => 'C'},
                {increment => +1, atom => 'N'},
            ],
            angle => 117.5,
        },
        {
            atoms => [
                {increment =>  0, atom => 'C'},
                {increment =>  0, atom => 'O'},
            ],
            angle => 121.5,
        },
    ],
    N => [
        {
            atoms => [
                {increment => 0, atom => 'CA'},
                {increment => 0, atom => 'C'},
            ],
            angle => 111.6
        },
    ],
    O => [
        {
            atoms => [
                {increment =>  0, atom => 'C'},
                {increment => +1, atom => 'N'},
            ],
            angle => 123.0
        },
    ],
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
   {increment => 0, atom => 'N'},
   {increment => 0, atom => 'CA'},
   {increment => 0, atom => 'C'},
   {increment => 1, atom => 'N'},
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
