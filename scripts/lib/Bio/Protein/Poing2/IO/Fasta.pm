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
    my @residues = @{Bio::Protein::Poing2::IO::Fasta::read_fasta('query.fasta')};

=head1 DESCRIPTION

This module provides a function for reading in a FASTA-formatted file and
returning a list of residues;

=head1 FUNCTIONS

=over

=item C<read_fasta($file)>: Read sequence from file C<$file>. C<$file> may
either be a filehandle or a file name.

=cut

sub read_fasta {
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

    my $res_num = 1;
    my %residues = ();
    while(my $ln = <$fh>){
        next if $ln =~ /^>/;
        chomp $ln;
        for(split //, $ln){
            $residues{$res_num++} =
                Bio::Protein::Poing2::Residue->new(type => $_)
        }
    }

    #Close filehandle if we opened it
    close $fh if ref $file ne 'GLOB';
    return \%residues;
}

=back

=cut
1;
