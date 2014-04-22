#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Temp;
use FindBin qw($Bin);
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

        #Pad rebuilt out to 80 chars so mkdssp doesn't crash
        open my $rebuilt_in, q{<}, "$tmp.rebuilt.pdb";
        open my $padded_out, q{>}, "$tmp.padded.pdb";
        while(my $ln = <$rebuilt_in>){
            chomp $ln;
            $ln .= " " x (80 - length $ln);
            print {$padded_out} $ln, "\n";
        }
        close $padded_out;
        close $rebuilt_in;

        system "mkdssp $tmp.padded.pdb > $tmp.dssp";
        system "$Bin/dssp2pdb $tmp.dssp $tmp.rebuilt.pdb > $tmp.new.pdb";

        open my $in, q{<}, "$tmp.new.pdb";
        my @rebuilt = grep {/^(ATOM|HELIX|SHEET)/} <$in>;
        close $in;
        unlink "$tmp.dssp", "$tmp.padded.pdb",
               "$tmp.rebuilt.pdb", "$tmp.new.pdb";
        @rebuilt;
    };

    printf "MODEL     % d\n", $model++;
    print @rebuilt;
    print "ENDMDL\n";

}
