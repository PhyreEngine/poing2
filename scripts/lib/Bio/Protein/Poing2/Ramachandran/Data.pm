package Bio::Protein::Poing2::Ramachandran::Data;
use strict;
use warnings;
use utf8;
use Moose;

=head1 NAME

Bio::Protein::Poing2::Ramachandran::Data - Find config data for Ramachandran constraints

=head1 SYNOPSIS

    #Guess config directory
    my $rama = Bio::Protein::Poing2::Ramachandran::Data->new();
    print $rama->string_repr;

    #Alternatively:
    my $rama = Bio::Protein::Poing2::Ramachandran::Data->new(
        dir => '/path/to/config/dir',
    );

=head1 DESCRIPTION

This class finds the boundary data for Ramachandran constraints. If the C<dir>
attribute is set, that directory is examined for files. Otherwise, first the
directory mentioned in the environment variable C<$RAMA_DATA> is examined (if
it is set), then C<$HOME/.poing2/rama/>. If none are set, die.

=cut

has dir => (
    is      => 'ro',
    default => undef,
    isa     => 'Maybe[Str]',
);

has files => (
    is       => 'ro',
    init_arg => undef,
    lazy     => 1,
    builder  => '_find_files',
);

use overload q{""} => \&string_repr;

sub _find_files {
    my ($self) = @_;

    my $files;
    if($self->dir && ($files = $self->_files_in_dir($self->dir))){
        return $files;
    }
    if($ENV{RAMA_DATA} && ($files = $self->_files_in_dir($ENV{RAMA_DATA}))){
        return $files;
    }
    if($files = $self->_files_in_dir('~/.phyrestorm/rama')){
        return $files;
    }
    die "Couldn't find Ramachandran data files.";
}

sub _files_in_dir {
    my ($self, $dir) = @_;

    #Expand path
    $dir = glob($dir);

    my %files = (
        ALANINE     => "$dir/boundary-ala-nosec.data",
        GENERAL     => "$dir/boundary-general-nosec.data",
        GLYCINE     => "$dir/boundary-gly-sym-nosec.data",
        PRE_PROLINE => "$dir/boundary-prepro.data",
        PROLINE     => "$dir/boundary-pro.data",
        ALPHA       => "$dir/boundary-alpha.data",
        BETA        => "$dir/boundary-beta.data",
    );
    for(values %files){
        return undef if !(-e $_);
    }
    return \%files;
}

sub string_repr {
    my ($self) = @_;

    my @lines = ();
    while(my ($type, $file) = each %{$self->files}){
        push @lines, "$type = $file\n";
    }
    return @lines;
}

sub TO_JSON {
    my ($self) = @_;
    return $self->files;
}

__PACKAGE__->meta->make_immutable;
1;
