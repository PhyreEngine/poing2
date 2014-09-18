#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 7;
use Test::Exception;

BEGIN {use_ok 'Bio::Protein::Poing2::Class'}

{
    package Foo;
    use Bio::Protein::Poing2::Class;
    use parent 'Bio::Protein::Poing2::Class';

    has mand => (is => 'ro', required => 1);
    has opt  => (is => 'rw', default  => sub {[0, 1, 2]});
}


throws_ok { Foo->new() }
    qr/^Mandatory argument `mand' not supplied/m,
    "Requires mand";
throws_ok { Foo->new(mand => 1, blah => 1) }
    qr/^Unknown attribute `blah'/m,
    "Rejects unknown attributes";

my $foo = Foo->new(mand => 'abc');
is($foo->mand, 'abc', 'Getter works');
throws_ok { $foo->mand('def') }
    qr/^Attribute `mand' is read-only/,
    'Rejects writing to read-only attributes';

is_deeply($foo->opt, [0, 1, 2], 'Default coderefs correct');
$foo->opt([1,2,3]);
is_deeply($foo->opt, [1, 2, 3], 'Setting rw attribute works');
