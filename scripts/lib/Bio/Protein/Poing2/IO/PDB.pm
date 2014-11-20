package Bio::Protein::Poing2::IO::PDB;
use strict;
use warnings;
use utf8;
use autodie;
use Bio::Protein::Poing2::Residue;
use Bio::Protein::Poing2::Atom;
use Bio::Protein::Poing2::Vector;

=head1 NAME

Bio::Protein::Poing2::IO::PDB - Read PDB files

=head1 SYNOPSIS

    use Bio::Protein::Poing2::IO::PDB;
    my %residues = %{ Bio::Protein::Poing2::IO::PDB::read_pdb('template.pdb')};

=head1 DESCRIPTION

C<Bio::Protein::Poing2::IO::PDB> provides the function C<read_pdb> for reading
a PDB file into a list of residues. Residues numbering is preserved as the keys
of the returned hashref.

=head1 FUNCTIONS

=over

=item C<read_pdb($file)>: Read a PDB file from the file C<$file>. C<$file> may
be a filehandle or a file name.

=back

=cut

sub read_pdb {
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
        next unless $ln =~ /^ATOM/;
        my $atom_name = substr $ln, 12, 4;
        my $res_num   = substr $ln, 22, 4;
        my $res_type  = substr $ln, 17, 3;
        my $x         = substr $ln, 30, 8;
        my $y         = substr $ln, 38, 8;
        my $z         = substr $ln, 46, 8;
        $res_num   =~ s/ //g;
        $atom_name =~ s/ //g;

        $residues{$res_num} ||= Bio::Protein::Poing2::Residue->new(
            index => $res_num,
            type  => $res_type,
        );
        my $atom = Bio::Protein::Poing2::Atom->new(
            name   => $atom_name,
            coords => Bio::Protein::Poing2::Vector->new(coords => [$x, $y, $z]),
            residue => $residues{$res_num},
        );
        push @{$residues{$res_num}->atoms}, $atom;
    }
    close $fh if ref $file ne 'GLOB';

    return \%residues;
}

1;
