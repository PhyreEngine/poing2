#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use autodie;
use constant PI => 4 * atan2(1, 1);

=head1 NAME

=head1 USAGE

pdb2springs.pl [PDB] [PDB2]

Generate a list of springs from the specified PDB files. If no PDB file is
specified, then standard input is read.

=head1 OPTIONS AND ARGUMENTS

=over

=item B<-h>, B<--help>

Display this help message.

=item B<--all-bb>

Include all backbone atoms instead of just a CA-sphere model.

=item B<-m>, B<--min-seq-sep> I<N>

Only place springs between residues at least I<N> residues apart. (Default: 3)

=item B<-M>, B<--max-seq-seqp> I<N>

Only place springs between residues at most I<N> residues apart. (Default: 20)

=item B<--overconstrain> I<N>

Place a maximum of I<N> springs between each residue pair. (Default: no limit)

=item B<--overconstrain-angle> I<N>

Place a maximum of I<N> torsion springs constraining each 4 residues. (Default:
no limit).

=item B<--between> I<ATOM>

Place springs between atom I<ATOM> of each residues. (Default: CA).

=item B<--query> I<FASTA>

FASTA file with the query.

=item B<--synth-time> I<ST>

Take I<ST> time units between synthesising each residue. (Default: 10)

=item B<-t>, B<--timestep> I<DT>

Each timestep is I<DT> time units. (Default: 0.1)

=item B<--const> I<K>

Use a spring constant of I<K>. (Default: 0.01)

=item B<--torsion-const> I<Kt>

Use a torsion spring constant of I<Kt>. (Default: 0.01)

=item B<--torsion-sep> I<t>

Form torsion springs betwee residues I<i>, I<i+t>, I<i+2t>, I<I+3t>. (Default:
1)

=item B<--cutoff> I<D>

Only allow springs to operate when they are less than I<D> Angstroms from their
equilibrium position. (Default: 10)

=item B<--max-dist> I<D>

Only place springs between atoms less than I<D> A apart. (Default: 20)

=item B<--no-torsion>

Disable torsion springs.

=item B<--no-linear>

Disable linear springs.

=item B<--no-sterics>

Disable steric effects.

=item B<--no-water>

Disable water battering.

=item B<--no-preamble>

Disable writing the preamble.

=item B<--fix>

Fix residues after the next residue is synthesised.

=item B<--drag> I<Cd>

Set drag coefficient to I<Cd>. This should be I<negative>. (Default: -0.1)

=back

=cut

our $atom_rec = "ATOM  % 5d %s %s A% 4d    %8.3f%8.3f%8.3f\n";
my %one2three = (
A=>'ALA', C=>'CYS', D=>'ASP', E=>'GLU', F=>'PHE', G=>'GLY', H=>'HIS', I=>'ILE',
K=>'LYS', L=>'LEU', M=>'MET', N=>'ASN', P=>'PRO', Q=>'GLN', R=>'ARG', S=>'SER',
T=>'THR', V=>'VAL', W=>'TRP', Y=>'TYR', Z=>'GLX', B=>'ASX',
);

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


my %options = (
    overconstrain         => undef,
    'overconstrain-angle' => undef,
    'min-seq-sep'         => 0,
    'max-seq-sep'         => 20,
    'between'             => 'CA',
    'synth-time'          => 10,
    'timestep'            => 0.1,
    'const'               => 0.01,
    'torsion-const'       => 0.01,
    'torsion-sep'         => 1,
    'cutoff'              => 10,
    'max-dist'            => 20,
    'drag'                => -0.1,
);
Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
    'all-bb',
    'min-seq-sep|m=i',
    'max-seq-sep|M=i',
    'no-torsion',
    'no-linear',
    'overconstrain=i',
    'overconstrain-angle=i',
    'between=s',
    'query|q=s',
    'synth-time=f',
    'timestep|t=f',
    'const|k=f',
    'torsion-const=f',
    'torsion-sep=i',
    'cutoff|c=f',
    'max-dist=f',
    'no-sterics',
    'no-water',
    'no-preamble',
    'add-bb',
    'fix',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

my @pdb_names = ();
my %pdbs = ();

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

my %residues = ();
while(<>){
    parse_atom(\%residues, $_) if /^ATOM/;
    if(eof(ARGV)){
        push @pdb_names, $ARGV;
        $pdbs{$ARGV} = {%residues};
        %residues = ();
    }
}

my $pairs = build_pairs(\@pdb_names, \%pdbs);
my $torsion = build_torsion_springs(\@pdb_names, \%pdbs);

if(not defined $options{'no-preamble'}){
    print "timestep = $options{timestep}\n"       if defined $options{timestep};
    print "synth_time = $options{'synth-time'}\n" if defined $options{'synth-time'};
    print "use_sterics = true\n"                  if !defined $options{'no-sterics'};
    print "use_water = true\n"                    if !defined $options{'no-sterics'};
    print "drag_coefficient = $options{'drag'}\n" if defined $options{'drag'};
    print "fix = true\n"                          if defined $options{'fix'};
    if(defined $options{'synth-time'} and defined $options{'query'}){
        printf "until = %d\n",
            (query_len($options{query}) + 5) * $options{'synth-time'};
    }
}

if(!$options{'no-pdb'} && $options{'query'}){
    if(!$options{'all-bb'}){
        print_pdb_ca_only($options{'query'});
    }else{
        print_pdb_all_bb($options{'query'});
    }
}

if(!$options{'no-linear'}){
    print "[Linear]\n";
    for my $pair_id(keys %{$pairs}){
        for my $spring(@{$pairs->{$pair_id}}){
            if($spring->{dist} < $options{'max-dist'}){
                printf "% 4d % 4s % 4d % 4s %8.6f %8.6f %8.6f\n",
                    $spring->{ra}{id}, $spring->{aa}{name},
                    $spring->{rb}{id}, $spring->{ab}{name},
                    $spring->{dist}, $options{const}, $options{cutoff};
            }
        }
    }
    #Print sidechain springs if possible
    if(!$options{'all-bb'} && $options{query}){
        open my $q_in, q{<}, $options{query};
        my $i = 1;
        while(<$q_in>){
            next if /^>/;
            chomp;
            for my $aa(split //){
                #Print backbone springs
                if($options{'add-bb'} && $i > 1){
                    printf "% 4d % 4s % 4d % 4s %8.6f %8.6f\n",
                        $i-1, 'CA', $i, 'CA', 3.8, 0.1;
                }
                next if $aa eq 'G';
                printf "% 4d % 4s % 4d % 4s %8.6f %8.6f\n",
                    $i, 'CA', $i, $one2three{$aa},
                    $bond_lens{$one2three{$aa}}, 0.1;

            }continue{$i++}
        }
        close $q_in;
    }

}

if(!$options{'no-torsion'}){
    print "[Torsion]\n";
    for my $id(sort {$a <=> $b} keys %{$torsion}){
        for my $spring(@{$torsion->{$id}}){
            printf "% 4d % 4d % 4d % 4d %8.5f %8.6f\n",
                $spring->{r1}{id},
                $spring->{r2}{id},
                $spring->{r3}{id},
                $spring->{r4}{id},
                $spring->{angle},
                $options{'torsion-const'};
        }
    }
}

sub print_pdb_ca_only {
    my ($query_file) = @_;
    print "[PDB]\n";
    open my $q_in, q{<}, $query_file;

    my $res = 0;
    while(<$q_in>){
        next if /^>/;
        chomp;

        for my $q(split //){
            printf $atom_rec,
                $res * 2 + 1, q{ CA }, $one2three{$q},
                $res + 1, 0, 0, 0;
            printf $atom_rec,
                $res * 2 + 2,
                qq{$one2three{$q} }, $one2three{$q},
                $res + 1, 0, 0, 0;
        }continue{$res++}
    }
    close $q_in;
}

sub print_pdb_all_bb {
    my ($query_file) = @_;
    print "[PDB]\n";
    open my $q_in, q{<}, $query_file;

    my @bb_atoms = (q{ N  }, q{ CA }, q{ C  }, q{ O  });
    my $res = 0;
    while(<$q_in>){
        next if /^>/;
        chomp;

        for my $q(split //){
            my $i = 1;
            for my $bb_atom(@bb_atoms){
                printf $atom_rec,
                    $res * 4 + $i, $bb_atom, $one2three{$q},
                    $res + 1, 0, 0, 0;
            }continue{ $i++ }
        }continue{$res++}
    }
    close $q_in;
}

sub print_query {
    my ($filename) = @_;

    print 'Sequence = ';
    open my $in, q{<}, $filename;
    while(<$in>){
        chomp;
        print unless /^>/;
    }
    close $in;
    print "\n";
}

sub build_pairs {
    my ($pdb_names, $pdbs) = @_;
    my @backbone = qw(N CA C O);

    my %pairs = ();
    for my $pdb_name(@{$pdb_names}){
        my $pdb = $pdbs->{$pdb_name};

        for my $i(sort {$a<=>$b} keys %{$pdb}){
            for my $j(sort {$a<=>$b} keys %{$pdb}){
                last if $j > $i;
                last if $i == $j && !$options{'all-bb'};
                #Only add pairs within the given bounds
                next if $i - $j > $options{'max-seq-sep'};
                next if $i - $j < $options{'min-seq-sep'};

                my $pair_id = "$i-$j";
                $pairs{$pair_id} ||= [];

                #Don't add another pair if we already have the required number
                #of constraints
                next if $options{overconstrain} &&
                        @{$pairs{$pair_id}} >= $options{overconstrain};

                if($options{'all-bb'}){
                    for my $k(0 .. $#backbone){
                        my $atom1 = $pdb->{$i}{atoms}{$backbone[$k]};
                        next unless $atom1;

                        for my $l(0 .. $#backbone){
                            next if $l >= $k && $i == $j;

                            my $atom2 = $pdb->{$j}{atoms}{$backbone[$l]};
                            next unless $atom2;

                            push @{$pairs{$pair_id}}, {
                                ra => $pdb->{$i},
                                rb => $pdb->{$j},
                                aa => $atom1,
                                ab => $atom2,
                                dist => dist($atom1, $atom2),
                            };
                        }
                    }
                }else{
                    my $atom1 = $pdb->{$i}{atoms}{$options{between}}
                        || $pairs{$pair_id}->{atoms}{'CA'};

                    my $atom2 = $pdb->{$j}{atoms}{$options{between}}
                        || $pairs{$pair_id}->{atoms}{'CA'};


                    push @{$pairs{$pair_id}}, {
                        ra => $pdb->{$i},
                        rb => $pdb->{$j},
                        aa => $atom1,
                        ab => $atom2,
                        dist => dist($atom1, $atom2),
                    };
                }
            }
        }
    }
    return \%pairs;
}

sub build_torsion_springs {
    my ($pdb_names, $pdbs) = @_;

    my $sep = $options{'torsion-sep'};
    my %springs = ();
    for my $pdb_name(@{$pdb_names}){
        my $pdb = $pdbs->{$pdb_name};

        for(my $i=1; $i<keys(%{$pdb}) - 2; $i++){
            next if !$pdb->{$i}        || !$pdb->{$i+1*$sep}
                 || !$pdb->{$i+2*$sep} || !$pdb->{$i+3*$sep};

            $springs{$i} ||= [];
            next if $options{'overconstrain-angle'} &&
                    @{$springs{$i}} >= $options{'overconstrain-angle'};

            my $spring = torsion_spring(
                $pdb->{$i},        $pdb->{$i+1*$sep},
                $pdb->{$i+2*$sep}, $pdb->{$i+3*$sep},
            );
            push @{$springs{$i}}, $spring;
        }
    }
    return \%springs;
}

sub dist {
    my ($a1, $a2) = @_;
    return sqrt(
        ($a1->{x} - $a2->{x})**2
        + ($a1->{y} - $a2->{y})**2
        + ($a1->{z} - $a2->{z})**2
    );
}

sub parse_atom {
    my ($residues, $line) = @_;

    my %atom = (
        id           => substr($line, 6, 5),
        name         => substr($line, 12, 4),
        conformation => substr($line, 16, 1),
        x            => substr($line, 30, 8),
        y            => substr($line, 38, 8),
        z            => substr($line, 46, 8),
    );
    s/ //g for values %atom;

    my %residue = (
        name     => substr($line, 17, 3),
        chain_id => substr($line, 21, 1),
        id       => substr($line, 22, 4),
        icode    => substr($line, 26, 1),
    );
    s/ //g for values %residue;

    $residues->{$residue{id}} ||= {
        %residue,
        atoms => {},
    };
    $residues->{$residue{id}}->{atoms}{$atom{name}} = \%atom;

}

sub torsion_spring {
    my ($r1, $r2, $r3, $r4) = @_;
    my $a1 = $r1->{atoms}{$options{between}};
    my $a2 = $r2->{atoms}{$options{between}};
    my $a3 = $r3->{atoms}{$options{between}};
    my $a4 = $r4->{atoms}{$options{between}};
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
    return {
        r1    => $r1,
        r2    => $r2,
        r3    => $r3,
        r4    => $r4,
        angle => $angle
    };
}

sub displacement {
    my ($a1, $a2) = @_;
    return {
        x => $a1->{x} - $a2->{x},
        y => $a1->{y} - $a2->{y},
        z => $a1->{z} - $a2->{z},
    };
}

sub cross {
    my ($a1, $a2) = @_;
    return {
        x => $a1->{y}*$a2->{z} - $a1->{z}*$a2->{y},
        y => $a1->{z}*$a2->{x} - $a1->{x}*$a2->{z},
        z => $a1->{x}*$a2->{y} - $a1->{y}*$a2->{x},
    };
}

sub dot {
    my ($a1, $a2) = @_;
    return $a1->{x}*$a2->{x} + $a1->{y}*$a2->{y} + $a1->{z}*$a2->{z};
}

sub mag {
    my ($a1) = @_;
    return sqrt(dot($a1, $a1));
}

sub div {
    my ($a, $s) = @_;
    return {
        x => $a->{x} / $s,
        y => $a->{y} / $s,
        z => $a->{z} / $s,
    };
}

sub query_len {
    my ($query) = @_;
    my $len = 0;
    open my $in, q{<}, $query;
    while(<$in>){
        next if /^>/;
        chomp;
        $len += length;
    }
    close $in;
    return $len;
}
