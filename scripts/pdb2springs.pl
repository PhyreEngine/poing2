#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use autodie;
use constant PI => 4 * atan2(1, 1);

=head1 NAME

pdb2springs.pl - Convert a bunch of PDB files to springs

=head1 USAGE

pdb2ssprings.pl B<FASTA> [B<PDB>] [B<PDB>...]

=cut

our $atom_rec = "ATOM  % 5d  %s %s A% 4d    %8.3f%8.3f%8.3f\n";
my %one2three = (
A=>'ALA', C=>'CYS', D=>'ASP', E=>'GLU', F=>'PHE', G=>'GLY', H=>'HIS', I=>'ILE',
K=>'LYS', L=>'LEU', M=>'MET', N=>'ASN', P=>'PRO', Q=>'GLN', R=>'ARG', S=>'SER',
T=>'THR', V=>'VAL', W=>'TRP', Y=>'TYR', Z=>'GLX', B=>'ASX',
);
my %three2one = ();
@three2one{values %one2three} = keys %one2three;

my %bond_lens = (
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

my %bb_len = (
    C  => {N => 1.3298236446806, O => 1.23203299355537},
    N  => {CA => 1.45999906298163},
    CA => {C => 1.53135743668343, CA => 3.8},
);

my %rama_phi_angles = (
    E => -90,
    H => -90,
    C =>  45,
);
my %rama_psi_angles = (
    E => 170,
    H => -20,
    C => 45,
);

my %rama_omega_angles = (
    E => 180,
    H => 180,
    C => 180,
);

my %options = (
    'min-seq-sep'     => 2,
    'max-seq-sep'     => 100,
    'bb-atom'         => 'CA',
    'timestep'        => 0.1,
    'synth-time'      => 10,
    'drag'            => -0.1,
    'post-time'       => 5,
    'const'           => 0.01,
    'cutoff'          => 10,
    'torsion-const'   => 0.001,
    'sidechain-const' => 0.1,
    'bb-const'        => 0.1,
);
Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
    'min-seq-sep=i',
    'max-seq-sep=i',
    'bb-atom=s',
    'all-bb',
    'timestep|t=f',
    'synth-time|T=f',
    'no-sterics',
    'no-water',
    'drag=f',
    'until|u=f',
    'timestep|d=f',
    'fix',
    'post-time=f',
    'no-sidechains',
    'const=f',
    'cutoff=f',
    'torsion-const=f',
    'no-preamble',
    'no-query',
    'no-linear',
    'no-torsion',
    'max-dist=f',
    'no-add-bb',
    'no-sc',
    'ss=s',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

my $backbone = ($options{'all-bb'})
    ? {CA => 1, C => 1, O => 1, N => 1}
    : {CA => 1};

my $backbone_linkage_prev = ($options{'all-bb'})
    ? {C => ['N' ]}
    : {$options{'bb-atom'} => [$options{'bb-atom'}]};

my $backbone_linkage_cur  = ($options{'all-bb'})
    ? {N => ['CA'], CA => ['C'], C => ['O']}
    : {};

my @phi = (
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
   {increment =>  0, atom => 'CA'},
   {increment =>  0, atom => 'C'},
);

my @psi = (
   {increment => -1, atom => 'N'},
   {increment => -1, atom => 'CA'},
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
);

my @omega = (
   {increment => -1, atom => 'CA'},
   {increment => -1, atom => 'C'},
   {increment =>  0, atom => 'N'},
   {increment =>  0, atom => 'CA'},
);

my @ca_ca_torsion = (
    {increment => -3, atom => 'CA'},
    {increment => -2, atom => 'CA'},
    {increment => -1, atom => 'CA'},
    {increment =>  0, atom => 'CA'},
);
my @ca_cb_torsion = (
    {increment => -1, atom => 'CA'},
    {increment => -1, atom => 'CB'},
    {increment =>  0, atom => 'CA'},
    {increment =>  0, atom => 'CB'},
);

my $fasta = shift || pod2usage('No fasta file supplied.');
my $query = load_query($fasta);
my $ss    = load_ss($options{ss}) if $options{ss};
my @pdbs  = @ARGV;
my %pdb_atoms = map {$_ => get_atoms($_)} @pdbs;
my $pairs = build_pairs(\%pdb_atoms);

#Set final time to a reasonable value if supplied
$options{until} ||= (@{$query} + $options{'post-time'})*$options{'synth-time'};

$pairs = filter_by_backbone($pairs, $backbone);


$pairs = filter_by_seq_sep($pairs,
    $options{'min-seq-sep'}, $options{'max-seq-sep'});

$pairs = filter_by_dist($pairs, $options{'max-dist'}) if $options{'max-dist'};

$pairs = add_sidechains($pairs, $query) unless $options{'all-bb'} || $options{'no-sc'};

$pairs = add_bb_springs($pairs, scalar(@{$query}),
    $backbone_linkage_prev, $backbone_linkage_cur) unless $options{'no-add-bb'};


#Begin printing output
print_preamble(\%options)                   unless $options{'no-preamble'};
print_query($query, $backbone, \%options)   unless $options{'no-query'};
print_linear_springs($pairs, \%options)     unless $options{'no-linear'};

if(!$options{'no-torsion'}){
    if($options{'all-bb'}){
        my $phi = build_dihedral_sets(\%pdb_atoms, scalar(@{$query}), \@phi);
        my $psi = build_dihedral_sets(\%pdb_atoms, scalar(@{$query}), \@psi);
        my $omega = build_dihedral_sets(\%pdb_atoms, scalar(@{$query}), \@omega);

        if($options{'ss'}){
            $phi = add_missing_dihedrals(
                $phi, scalar(@{$query}), \@phi, $ss, \%rama_phi_angles);
            $psi = add_missing_dihedrals(
                $psi, scalar(@{$query}), \@psi, $ss, \%rama_psi_angles);
            $omega = add_missing_dihedrals(
                $omega, scalar(@{$query}), \@omega, $ss, \%rama_omega_angles);
        }

        my $fourmers = [@{$phi}, @{$psi}, @{$omega}];
        print_torsion_springs($fourmers, \%options);
    }else{
        my $torsion = build_dihedral_sets(\%pdb_atoms, scalar(@{$query}), \@ca_ca_torsion);

        #Build angles for sidechains. Do that by working with CB atoms and then
        #rename them appropriately.
        my $sc = build_dihedral_sets(\%pdb_atoms, scalar(@{$query}), \@ca_cb_torsion);
        for my $fourmer(@{$sc}){
            my $cb1 = $fourmer->{fourmer}[1];
            my $cb3 = $fourmer->{fourmer}[3];
            $cb1->{atom} = $one2three{$query->[$cb1->{res} - 1]};
            $cb3->{atom} = $one2three{$query->[$cb3->{res} - 1]};
        }
        print_torsion_springs([@{$torsion}, @{$sc}], \%options);
    }
}

=begin comment

 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

=end comment

=cut

#Print the preamble
sub print_preamble {
    my ($options) = @_;

    my %opt_map = (
        timestep         => 'timestep',
        synth_time       => 'synth-time',
        use_sterics      => '!no-sterics',
        use_water        => '!no-sterics',
        drag_coefficient => 'drag',
        fix              => 'fix',
        until            => 'until',
    );

    my %boolean = map {$_ => 1} qw(no-sterics fix);
    while(my ($param, $opt_param) = each %opt_map){
        my ($op) = $opt_param =~ /^!?(.*)$/;

        my $bool = exists $options->{$op};
        #Invert if necessary
        $bool = !$bool if $opt_param =~ /^!/;

        #The value is the value in $options if it is not boolean. Otherwise it
        #must be "true" or "false".
        my $val;
        if($boolean{$op}){
            $val = $bool ? 'true' : 'false';
        }else{
            $val = $options{$op};
        }
        print "$param = $val\n" if $bool;
    }
}

#Add sidechain atom to each residue
sub add_sidechains {
    my ($pairs, $seq) = @_;

    my @sc = ();
    my $i = 1;
    for my $aa(@{$seq}){
        next if $aa eq 'G';
        push @sc, {
            i => $i,
            j => $i,
            atom_i => 'CA',
            atom_j => $one2three{$aa},
            dist => $bond_lens{$one2three{$aa}},
            const => $options{'sidechain-const'},
        };
    }continue{ $i++ }

    push @{$pairs}, @sc;
    return $pairs;
}

#Print query in PDB-ish format
sub print_query {
    #$query:    arrayref of single character AAs
    #$backbone: hashref with the keys as atom times
    my ($query, $backbone, $options) = @_;

    my $atoms_per_res = keys %{$backbone};
    $atoms_per_res++ unless $options{'all-bb'} || $options{'no-sidechains'};

    print "[PDB]\n";
    my $res = 0;
    for my $aa(@{$query}){
        my $i = 1;
        for my $atom(keys %{$backbone}){
            printf $atom_rec,
                $res * $atoms_per_res + $i,
                sprintf('%-3s', $atom),
                $one2three{$aa},
                $res + 1,
                0, 0, 0;
        }continue{ $i++ }

        #Print sidechain
        if($aa ne 'G' && (!$options{'all-bb'} && !$options{'no-sidechains'})){
            printf $atom_rec,
                $res * $atoms_per_res + $i,
                sprintf('%-3s', $one2three{$aa}),
                $one2three{$aa},
                $res + 1,
                0, 0, 0;
        }
    }continue{ $res++ }
}

#Print linear springs
sub print_linear_springs {
    #pairs:   Arrayref of pairs
    #options: Hashref of options
    my ($pairs, $options) = @_;
    print "[Linear]\n";
    for my $pair(@{$pairs}){
        printf "% 4d % 4s % 4d % 4s %8.6f %8.6f %8.6f\n",
            $pair->{i}, $pair->{atom_i},
            $pair->{j}, $pair->{atom_j},
            $pair->{dist},
            ($pair->{const} // $options->{const}),
            $options->{cutoff};
    }
}

#Print torsion springs
sub print_torsion_springs {
    my ($fourmers, $options) = @_;
    print "[Torsion]\n";
    for my $fourmer(@{$fourmers}){
        printf "% 4d % 4s % 4d % 4s % 4d % 4s % 4d % 4s %8.5f %8.6f\n",
            $fourmer->{fourmer}[0]{res}, $fourmer->{fourmer}[0]{atom},
            $fourmer->{fourmer}[1]{res}, $fourmer->{fourmer}[1]{atom},
            $fourmer->{fourmer}[2]{res}, $fourmer->{fourmer}[2]{atom},
            $fourmer->{fourmer}[3]{res}, $fourmer->{fourmer}[3]{atom},
            $fourmer->{angle},
            $options{'torsion-const'};
    }
}


#Add missing dihedral angles from Ramachandran space using the (possibly
#estimated) estimated secondary structure of each residue. This should be
#called immediately after build_dihedral_sets with the same link argument,
#because it will not check to see which atoms are actually in the fourmer.
#
#The first argument must be a list of fourmers.
#
#The second argument is the query length.
#
#The third argument should describe all atoms between which to add atoms,
#similarly to build_dihedral_sets.
#
#The fourth argument is an arrayref of secondary structure states.
#
#The fifth argument is a hashref of SS state mapping to dihedral angles, e.g.
#from a Ramachandran plot.

sub add_missing_dihedrals {
    my ($fourmers, $query_length, $link, $ss, $rama_angles) = @_;
    my %unconnected = ();

    #Map ss to residue numbers
    my %ss = ();
    @ss{1 .. $query_length} = @{$ss};

    #Find out where we need to start. That is, if we need to look back 3 spaces
    #we should start at residue 3.

    for my $i(2 .. $query_length){
        $unconnected{$i} = 1;
    }
    for my $fourmer(@{$fourmers}){
        my ($ca) = grep {$_->{atom} eq 'CA'} @{$fourmer->{fourmer}};
        delete $unconnected{$ca->{res}};
    }

    for my $res(keys %unconnected){
        my %dihedral = (angle => $rama_angles->{$ss{$res}});
        my @fourmer = map {
            {res => $res + $_->{increment}, atom => $_->{atom}}
        } @{$link};
        $dihedral{fourmer} = \@fourmer;
        push @{$fourmers}, \%dihedral;
    }
    return $fourmers;
}

#Scan through the atoms and get dihedral angles for each.
#
#The first argument should be a hashref of pdbs and arrayrefs of atoms, as
#returned by get_atoms.
#
#The second argument is the query length.
#
#The third argument describes which atoms should be linked. This must be an
#arrayref of hashrefs describing an increment to be added to the current
#residue number and the atom name. For example (phi):
#
# [
#   {increment => -1, atom => 'C'},
#   {increment =>  0, atom => 'N'},
#   {increment =>  0, atom => 'CA'},
#   {increment =>  0, atom => 'C'},
# ]
#
#The return value is an arrayref of hashrefs, each containing the keys "angle"
#and "fourmer". The "fourmer" key is an arrayref of the atoms in the fourmer.

sub build_dihedral_sets {
    my ($pdb_atoms, $query_length, $links) = @_;

    my @fourmers = ();
    while(my ($pdb, $atoms) = each %{$pdb_atoms}){
        #Generate list of atoms
        my %all_atoms = ();
        for my $atom(@{$atoms}){
            $all_atoms{$atom->{res}} ||= {};
            $all_atoms{$atom->{res}}->{$atom->{atom}} = $atom;
        }

        #Go over all the atoms and build the 4mers
        RESIDUE: for my $i(1 .. $query_length){
            my @in_fourmer = ();
            for my $desc(@{$links}){
                my $inc = $desc->{increment};
                my $atom = $all_atoms{$i + $inc}->{$desc->{atom}};
                next RESIDUE if !$atom;
                push @in_fourmer, $atom;
            }
            my $angle = dihedral_angle(map {$_->{coords}} @in_fourmer);
            push @fourmers, {fourmer => \@in_fourmer, angle => $angle};
        }
    }
    return \@fourmers;
}

#Get dihedral angle of four atoms
sub dihedral_angle {
    my ($a1, $a2, $a3, $a4) = @_;

    my $b1 = displacement($a2, $a1);
    my $b2 = displacement($a3, $a2);
    my $b3 = displacement($a4, $a3);

    my $cross_b1_b2 = cross($b1, $b2);
    my $cross_b2_b3 = cross($b2, $b3);
    my $cross1 = cross($cross_b1_b2, $cross_b2_b3);
    my $b2_mag = mag($b2);
    my $norm_b2 = div($b2, $b2_mag);
    my $y = dot($cross1, $norm_b2);
    my $x = dot($cross_b1_b2, $cross_b2_b3);

    my $angle = atan2($y, $x) * 180 / PI;
    return $angle;
}

sub displacement {
    my ($a1, $a2) = @_;
    return [
        $a1->[0] - $a2->[0],
        $a1->[1] - $a2->[1],
        $a1->[2] - $a2->[2],
    ];
}

sub cross {
    my ($a1, $a2) = @_;
    return [
        $a1->[1]*$a2->[2] - $a1->[2]*$a2->[1],
        $a1->[2]*$a2->[0] - $a1->[0]*$a2->[2],
        $a1->[0]*$a2->[1] - $a1->[1]*$a2->[0],
    ];
}

sub dot {
    my ($a1, $a2) = @_;
    return $a1->[0]*$a2->[0] + $a1->[1]*$a2->[1] + $a1->[2]*$a2->[2];
}

sub mag {
    my ($a1) = @_;
    return sqrt(dot($a1, $a1));
}

sub div {
    my ($a, $s) = @_;
    return [
        $a->[0] / $s,
        $a->[1] / $s,
        $a->[2] / $s,
    ];
}


#Add springs between adjacent backbone atoms if they are not actually
#connected.
#
#The first argument is the arrayref of pairs.
#
#The second argument is the query length.
#
#The third argument describes which atoms from i-1th residue will be linked to
#the ith residue; for example, {CA => ['CA']}.
#
#The fourth argument describes which atoms from the ith residue will be linked
#to the ith residue; for example, {N => ['CA'], CA => ['C'], C => ['O']}

sub add_bb_springs {
    my ($pairs, $query_len, $prev_con, $curr_con) = @_;

    #Begin by assuming all are unconnected
    my %unconnected = ();
    for my $i(1 .. $query_len){
        #Intra-residue
        for my $start_atom(keys %{$curr_con}){
            for my $end_atom(@{$curr_con->{$start_atom}}){
                $unconnected{"$i-$start_atom:$i-$end_atom"} = 1;
            }
        }

        #To previous residue
        next if $i == 1;
        for my $start_atom(keys %{$prev_con}){
            my $prev_res = $i-1;
            for my $end_atom(@{$prev_con->{$start_atom}}){
                $unconnected{"$prev_res-$start_atom:$i-$end_atom"} = 1;
            }
        }
    }

    for my $pair(@{$pairs}){
        my $ai = $pair->{atom_i};
        my $aj = $pair->{atom_j};
        delete $unconnected{"$pair->{i}-$ai:$pair->{j}-$aj"};
    }

    for my $con(keys %unconnected){
        my ($i, $atom_i, $j, $atom_j) = $con =~ /^(\d+)-(\w+):(\d+)-(\w+)$/;
        my $len = $bb_len{$atom_i}->{$atom_j};
        push @{$pairs}, {
            i => $i,
            j => $j,
            atom_i => $atom_i,
            atom_j => $atom_j,
            dist   => $len,
            const  => $options{'bb-const'},
        };
    }
    return $pairs;
}

#Filter by sequence separation
sub filter_by_seq_sep {
    my ($pairs, $min, $max) = @_;
    my @new_pairs = ();
    for my $pair(@{$pairs}){
        next if abs($pair->{i} - $pair->{j}) < $min;
        next if abs($pair->{i} - $pair->{j}) > $max;
        push @new_pairs, $pair;
    }
    return \@new_pairs;
}

#Filter by maximum distance
sub filter_by_dist {
    #$pairs: Arrayref of pairs
    #$bb:    Hashref of atoms
    my ($pairs, $max) = @_;
    return [grep {$_->{dist} <= $max} @{$pairs}];
}

#Filter the list of pairs by the backbone atoms required.
sub filter_by_backbone {
    #$pairs: Arrayref of pairs
    #$bb:    Hashref of atoms
    my ($pairs, $bb) = @_;
    my @new_pairs = ();
    for my $pair(@{$pairs}){
        next if !$bb->{$pair->{atom_i}} || !$bb->{$pair->{atom_j}};
        push @new_pairs, $pair;
    }
    return \@new_pairs;
}

#Read a query sequence from a FASTA file. This function returns an arrayref of
#letters.
sub load_query {
    my ($fasta) = @_;

    my @seq = ();
    open my $in, q{<}, $fasta;
    while(<$in>){
        next if /^>/;
        chomp;
        push @seq, split //;
    }
    close $in;
    return \@seq;
}

#Load a psipred ss2-style secondary structure file and return an arrayref of
#single-letter codes representing the secondary structure at each residue. It
#is assumed that all residues appear in the file and that all residues appear
#in order.
sub load_ss {
    my ($ss) = @_;
    my @secondary_structure = ();
    open my $in, q{<}, $ss;
    while(my $ln = <$in>){
        chomp $ln;
        next if $ln =~ /^#/;
        next if $ln =~ /^\s*$/;
        my $state = substr $ln, 7, 1;
        push @secondary_structure, $state;
    }
    close $in;
    return \@secondary_structure;
}


#Get all pairs from all PDB files. The argument should be a hashref of
#{$pdb_name => $atoms_arrayref}.
#
#The return value looks like:
#
# [
#   {
#       i => 1,         j => 2,
#       atom_i => 'CA', atom_j => 'CA',
#       dist => 3.0,    pdb => 'foo.pdb'
#   },
#   ...
# ]
sub build_pairs {
    my ($pdb_atoms) = @_;

    my @pairs = ();

    while(my ($pdb, $atoms) = each %{$pdb_atoms}){
        for my $a_i(@{$atoms}){
            for my $a_j(@{$atoms}){
                next if $a_j >= $a_i;
                push @pairs, {
                    i      => $a_i->{res},
                    j      => $a_j->{res},
                    atom_i => $a_i->{atom},
                    atom_j => $a_j->{atom},
                    dist   => dist($a_i->{coords}, $a_j->{coords}),
                    pdb    => $pdb,
                };
            }
        }
    }
    return \@pairs;
}

#Get atoms from PDB file
sub get_atoms {
    my ($pdb) = @_;
    my @atoms = ();

    open my $in, q{<}, $pdb;
    while(my $ln = <$in>){
        next unless $ln =~ /^ATOM/;
        my $atom = substr $ln, 12, 4;
        my $res  = substr $ln, 22, 4;
        my $x    = substr $ln, 30, 8;
        my $y    = substr $ln, 38, 8;
        my $z    = substr $ln, 46, 8;
        $atom =~ s/ //g;
        push @atoms, {res => int($res), atom => $atom, coords => [$x, $y, $z]};
    }
    close $in;
    return \@atoms;
}

#Get the list of pairs from a PDB. This includes all pairs; filtering should be
#applied later;

#Get distance between two coords
sub dist {
    my ($a, $b) = @_;
    return sqrt(
        ($a->[0] - $b->[0]) ** 2 +
        ($a->[1] - $b->[1]) ** 2 +
        ($a->[2] - $b->[2]) ** 2
    );
}

sub print_rama_plot {
    my ($pdb_atoms, $query) = @_;
    my $psi = build_dihedral_sets($pdb_atoms, scalar(@{$query}), \@psi);
    my $phi = build_dihedral_sets($pdb_atoms, scalar(@{$query}), \@phi);

    my %phi_psi_map = ();
    for my $fourmer(@{$psi}){
        my ($ca) = grep {$_->{atom} eq 'CA'} @{$fourmer->{fourmer}};
        $phi_psi_map{$ca->{res}} ||= {};
        $phi_psi_map{$ca->{res}}->{psi} = $fourmer->{angle};
    }
    for my $fourmer(@{$phi}){
        my ($ca) = grep {$_->{atom} eq 'CA'} @{$fourmer->{fourmer}};
        $phi_psi_map{$ca->{res}} ||= {};
        $phi_psi_map{$ca->{res}}->{phi} = $fourmer->{angle};
    }
    print "phi\tpsi\n";
    for(values %phi_psi_map){
        next unless $_->{phi} && $_->{psi};
        print "$_->{phi}\t$_->{psi}\n";
    }
}

