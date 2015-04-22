package Bio::Protein::Poing2::IO::Psipred;
use strict;
use warnings;
use utf8;
use autodie;
use Bio::Protein::Poing2::Residue;

=head1 NAME

Bio::Protein::Poing2::IO::Psipred - Read PSIPRED ss2-style file

=head1 SYNOPSIS

    use Bio::Protein::Poing2::IO::Psipred;
    my @residues = @{Bio::Protein::Poing2::IO::Fasta::read_ss2('query.ss2')};

=head1 DESCRIPTION

This module provides a function for reading in an SS2-formatted file and
returning a hash of residues;

=head1 FUNCTIONS

=over

=item C<read_ss2($file)>: Read secondarty structure  from file C<$file>.
C<$file> may either be a filehandle or a file name.

=cut

sub read_psipred {
    my ($file) = @_;
    my $fh = undef;
    if(ref $file eq 'GLOB'){
        $fh = $file;
    }else{
        #We're not going to forget to destroy this
        ##no critic (InputOutput::RequireBriefOpen)
        open $fh, q{<}, $file;
        ##use critic
    }

    my %residues = ();
    while(my $ln = <$fh>){
        chomp $ln;
        next if $ln =~ /^#/;
        next if $ln =~ /^\s*$/;
        my $index    = substr $ln, 0, 4;
        my $state    = substr $ln, 7, 1;
        my $res_type = substr $ln, 5, 1;
        $index =~ s/ //g;

        my $res = Bio::Protein::Poing2::Residue->new(
            type     => $res_type,
            ss_state => $state,
            index    => $index,
        );
        $residues{$index} = $res;
    }

    #Close filehandle if we opened it
    close $fh if ref $file ne 'GLOB';
    return \%residues;
}

=back

=cut
1
