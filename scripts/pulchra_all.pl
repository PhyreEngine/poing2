#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Temp;
use autodie;

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

my $model = 1;
my $tmp = File::Temp->new;
$/ = "ENDMDL\n";
while(<>){
    $tmp->truncate(0);
    $tmp->flush;
    my @atoms = grep {/^ATOM.* CA /} split "\n";
    next unless @atoms > 2;

    print {$tmp} join("\n", @atoms);
    $tmp->flush;

    system "pulchra $tmp >/dev/null";
    next if $? != 0;
    my @rebuilt = do {
        local $/ = "\n";
        open my $in, q{<}, "$tmp.rebuilt.pdb";
        my @rebuilt = grep {/^ATOM/} <$in>;
        close $in;
        @rebuilt;
    };

    printf "MODEL     % d\n", $model++;
    print @rebuilt;
    print "ENDMDL\n";

}
