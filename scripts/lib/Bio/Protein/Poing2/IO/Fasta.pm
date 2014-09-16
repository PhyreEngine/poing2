package Bio::Protein::Poing2::IO::Fasta;
use strict;
use warnings;
use utf8;
use autodie;
use Bio::Protein::Poing2::Residue;

=head1 NAME

Bio::Protein::Poing2::IO::Fasta - Read FASTA files

=head1 SYNOPSIS

    use Bio::Protein::Poing2::IO::Fasta;
    my @residues = @{ Bio::Protein::Poing2::IO::Fasta::read('query.fasta') };

=head1 DESCRIPTION

This module provides a function for reading in a FASTA-formatted file and
returning a list of residues;

=head1 FUNCTIONS

=over

=item C<read($file)>: Read sequence from file C<$file>. C<$file> may either be
a filehandle or a file name.

=cut

sub read {
    my ($file) = @_;
    my $fh = undef;
    if(ref $file eq 'GLOB'){
        $fh = $file;
    }else{
        open $fh, q{<}, $file;
    }

    my @residues = ();
    while(my $ln = <$fh>){
        next if $ln =~ /^>/;
        chomp $ln;
        push @residues, map {Bio::Protein::Poing2::Residue->new($_)}
                        split //, $ln;
    }

    #Close filehandle if we opened it
    close $fh if ref $file ne 'GLOB';
    return \@residues;
}

=back

=cut
1;
