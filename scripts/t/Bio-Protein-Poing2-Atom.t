#!/usr/bin/env perl
use strict;
use warnings;
use Math::Vector::Real;
use Test::More tests => 7;
use Test::Exception;

BEGIN {use_ok 'Bio::Protein::Poing2::Atom'}

my $atom = Bio::Protein::Poing2::Atom->new(
    name  => 'CA',
    index => 1,
);
is($atom->name, 'CA', 'Constructor sets name');
is($atom->x,    0,    'x coord defaults to zero');
is($atom->y,    0,    'y coord defaults to zero');
is($atom->z,    0,    'z coord defaults to zero');

$atom->x(1);
$atom->y(2);
$atom->z(3);
ok($atom->coords == V(1, 2, 3), 'Individual setters work');
$atom->coords(V(3, 2, 1));
ok($atom->coords == V(3, 2, 1), 'Group setter works');
