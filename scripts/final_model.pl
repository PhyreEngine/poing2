#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

=head1 NAME

=head1 USAGE

=head1 OPTIONS AND ARGUMENTS

=cut

my %options = ();
Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

$/ = "ENDMDL\n";

my $last = "";
while(<>){
    $last = $_;
}
print $last;
