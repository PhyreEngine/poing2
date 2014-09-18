#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 4;

BEGIN {use_ok 'Bio::Protein::Poing2::IO::PDB'}

my $pdb = <<EOPDB;
ATOM      1  N   LYS A   5      25.181  39.795  -5.151  1.00 30.25           N  
ATOM      2  CA  LYS A   5      24.955  38.412  -5.598  1.00 23.34           C  
ATOM      3  C   LYS A   5      24.154  38.507  -6.878  1.00 18.80           C  
ATOM      4  O   LYS A   5      23.741  39.593  -7.299  1.00 22.27           O  
ATOM      5  CB  LYS A   5      24.175  37.667  -4.526  1.00 29.50           C  
ATOM      6  CG  LYS A   5      22.827  38.272  -4.201  1.00 39.44           C  
ATOM      7  CD  LYS A   5      22.358  37.862  -2.793  1.00 51.81           C  
ATOM      8  CE  LYS A   5      20.836  37.777  -2.718  1.00 58.22           C  
ATOM      9  NZ  LYS A   5      20.318  36.822  -3.741  1.00 60.79           N  
ATOM     10  N   LYS A   6      23.911  37.371  -7.492  1.00 18.03           N  
ATOM     11  CA  LYS A   6      23.040  37.316  -8.632  1.00 18.43           C  
ATOM     12  C   LYS A   6      21.621  36.986  -8.134  1.00 19.54           C  
ATOM     13  O   LYS A   6      21.431  36.383  -7.096  1.00 18.32           O  
ATOM     14  CB  LYS A   6      23.559  36.280  -9.628  1.00 19.61           C  
ATOM     15  CG  LYS A   6      24.945  36.594 -10.192  1.00 19.22           C  
ATOM     16  CD  LYS A   6      25.374  35.511 -11.174  1.00 19.17           C  
ATOM     17  CE  LYS A   6      26.860  35.549 -11.572  1.00 20.41           C  
ATOM     18  NZ  LYS A   6      27.213  34.378 -12.422  1.00 22.29           N  
ATOM     19  N   PRO A   7      20.603  37.457  -8.851  1.00 23.13           N  
ATOM     20  CA  PRO A   7      19.223  37.177  -8.426  1.00 22.43           C  
ATOM     21  C   PRO A   7      18.781  35.746  -8.712  1.00 20.89           C  
ATOM     22  O   PRO A   7      19.381  35.009  -9.505  1.00 19.37           O  
ATOM     23  CB  PRO A   7      18.374  38.124  -9.237  1.00 27.06           C  
ATOM     24  CG  PRO A   7      19.156  38.262 -10.453  1.00 26.52           C  
ATOM     25  CD  PRO A   7      20.635  38.285 -10.047  1.00 25.70           C  
ATOM     26  N   GLY A   8      17.677  35.425  -8.070  1.00 18.76           N  
ATOM     27  CA  GLY A   8      17.073  34.128  -8.187  1.00 14.97           C  
ATOM     28  C   GLY A   8      17.462  33.206  -7.061  1.00 14.49           C  
ATOM     29  O   GLY A   8      18.317  33.503  -6.240  1.00 15.92           O  
ATOM     30  N   LEU A   9      16.766  32.073  -7.015  1.00 12.49           N
EOPDB

open my $pdb_in, q{<}, \$pdb;
my %res = %{Bio::Protein::Poing2::IO::PDB::read_pdb($pdb_in)};
close $pdb_in;

is(keys %res, 5, 'Read 5 residues');
is(@{$res{5}->atoms}, 9, '9 atoms in residue 5');
is(@{$res{9}->atoms}, 1, '1 atom in residue 9');
