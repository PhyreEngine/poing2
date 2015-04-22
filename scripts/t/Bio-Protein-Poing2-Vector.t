#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 9;
use Test::Exception;

BEGIN {*Poing2:: = *Bio::Protein::Poing2::}
BEGIN {use_ok 'Bio::Protein::Poing2::Vector'}

my $a = Poing2::Vector->new(coords => [1,2,3]);
my $b = Poing2::Vector->new(coords => [3,2,1]);
my $c = Poing2::Vector->new(coords => [3,2,1,0]);

my $r1 = Poing2::Vector->new(coords => [-2, 0,  2]);
my $r2 = Poing2::Vector->new(coords => [ 2, 0, -2]);

throws_ok { $a + $c } qr/^Dimensionality incorrect/, 'Requires equal dims';

ok($b - $a == Poing2::Vector->new(coords => [2, 0, -2]), 'Subtraction works');
ok($a + $b == Poing2::Vector->new(coords => [4, 4, 4]),  'Addition works');
ok($a * 3  == Poing2::Vector->new(coords => [3, 6, 9]),  'Multiplication works');
ok(($a * 3) / 3 == $a,                                   'Division works');
ok(-$a == $a * -1,                                       'Unary minus works');
ok($a . $b == 3 + 4 + 3,                                 'Dot product');
ok(
    Poing2::Vector->new(coords => [1, 0, 0])
    x
    Poing2::Vector->new(coords => [0, 1, 0])
    ==
    Poing2::Vector->new(coords => [0, 0, 1]),
    'Cross product works');
