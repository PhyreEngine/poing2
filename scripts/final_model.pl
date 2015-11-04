#!/usr/bin/env perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use autodie;

=head1 NAME

final_model.pl - Extract the final model from a multi-model PDB file

=head1 USAGE

B<final_model.pl> [B<PDB>]

If the argument B<PDB> is not given, standard input is read. Please note that
reading from standard input is much slower than reading from a file, as we
cannot seek to the end of a stream.

=cut

my %options = ();
Getopt::Long::Configure(qw(bundling no_ignore_case));
GetOptions(\%options,
    'help|h',
) or pod2usage(2);
pod2usage(-verbose => 2, -noperldoc => 1, -exitval => 1) if $options{help};

#If we're given a file, we can seek to the end and go backwards

if(@ARGV == 1){
    my $pdb = shift;
    open my $in, q{<}, $pdb;

    my $buffer = q{};

    #Seek to end of file
    seek $in, -1, 2;
    #Keep moving back 80 bytes (the length of an ATOM line) and prepend to buffer.
    while($buffer !~ /^MODEL/m){
        seek $in, -80, 1;
        my $tmp = q{};
        my $nread = read $in, $tmp, 80;
        seek $in, -$nread, 1;
        $buffer = $tmp . $buffer;
    }
    my @lines = split qq{\n}, $buffer;

    for(@lines){
        print $_, "\n" if /^MODEL/../^ENDMDL/;
    }
}else{
    local $/ = "ENDMDL\n";
    my $last = q{};
    $last = $_ while <>;
    print $last;
}
