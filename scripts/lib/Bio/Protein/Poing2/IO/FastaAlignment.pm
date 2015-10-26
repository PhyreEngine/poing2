package Bio::Protein::Poing2::IO::FastaAlignment;
use strict;
use warnings;
use utf8;
use autodie;
use Bio::Protein::Poing2::Residue;

=head1 NAME

Bio::Protein::Poing2::IO::FastaAlignment - Read FASTA alignment

=head1 SYNOPSIS

    use Bio::Protein::Poing2::IO::FastaAlignment;
    my %residues = %{Bio::Protein::Poing2::IO::FastaAlignment::read_fasta('aln.fasta')};

=head1 DESCRIPTION

This module provides a function for reading in a FASTA-formatted file and
returning a hash of residues. Consider the following alignment:

    >query
    AAG--AG
    >template
    TA----A

The result will be a hashref. Each key will be a residue number, and the value
will be a hashref containing a C<from> key and a C<to> key. The values of these
keys will be L<Bio::Protein::Poing2::Residue> objects of the correct type. The
above query would look like this (with the Residue objects represented as
single letters):

    {
        1 => {from => A, to => T},
        2 => {from => A, to => A},
        3 => {from => G, to => undef},
        4 => {from => A, to => undef},
        5 => {from => G, to => A},
    }

=head1 FUNCTIONS

=over

=item C<read_fasta($file)>: Read alignment from file C<$file>. C<$file> may
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

    my $seq_num = -1;
    my @sequences = (q{}, q{});

    #Build string sequences
    while(my $ln = <$fh>){
        if($ln =~ /^>/){
            $seq_num++;
            next;
        }
        chomp $ln;
        $sequences[$seq_num] .= $ln;
    }

    #Number residues by the query, starting from 1
    for(my $i=0; $i < length($sequences[0]); $i++){
        my $qseq = substr($sequences[0], $i, 1);
        my $tseq = substr($sequences[1], $i, 1);

        next if $qseq eq '-';
        my $from = Bio::Protein::Poing2::Residue->new(
            index => $i + 1,
            type  => $qseq,
        );
        my $to = undef;
        if($tseq ne '-'){
            $to = Bio::Protein::Poing2::Residue->new(
                index => $i + 1,
                type  => $tseq,
            );
        }
        $residues{$i+1} = {
            from => $from,
            to => $to,
        };
    }

    #Close filehandle if we opened it
    close $fh if ref $file ne 'GLOB';
    return \%residues;
}

=back

=cut
1;
