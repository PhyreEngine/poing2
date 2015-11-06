package Bio::Protein::Poing2;
use strict;
use warnings;
use Moose;
use Bio::Protein::Poing2::Ramachandran::Data;

=head1 Attributes

=over

=item C<query>

Query object.

=cut

has query => (is => 'ro', isa => 'Bio::Protein::Poing2::Query', required => 1);

=item C<templates>

Arrayref of template objects.

=cut

has templates => (is => 'ro', isa => 'ArrayRef[Bio::Protein::Poing2::Template]', required => 1);

=item C<spring_filters>

Filter linear springs with these filters.

=cut

has spring_filters => (
    is      => 'ro',
    isa     => 'ArrayRef[Bio::Protein::Poing2::Filter::Pair]',
    default => sub { [] },
);


=item C<rama_data> [optional]

Optional Ramachandran data object.

=cut

has rama_data => (
    is      => 'ro',
    isa     => 'Bio::Protein::Poing2::Ramachandran::Data',
    lazy    => 1,
    builder => '_build_rama_data',
);


use overload q{""} => \&to_str;

sub _build_rama_data {
    my ($self) = @_;
    return Bio::Protein::Poing2::Ramachandran::Data->new;
}

sub to_str {
    my ($self) = @_;

    my @lines;
    push @lines, "[PDB]\n";
    for my $i(sort {$a <=> $b} keys %{$self->query->backbone}){
        push @lines, $self->query->residues->{$i}->string_repr;
    }

    push @lines, "[Linear]\n";
    for my $pair(@{$self->query->backbone_springs}){
        push @lines, $pair->string_repr;
    }
    for my $template(@{$self->templates}){
        my $pairs = $template->pairs;
        $pairs = $_->filter($pairs) for @{$self->spring_filters};
        for my $pair(@{$pairs}){
            push @lines, $pair->string_repr;
        }
    }

    push @lines, "[Angle]\n";
    for my $angle(@{$self->query->angles}){
        push @lines, $angle->string_repr;
    }

    push @lines, "[Torsion]\n";
    for my $template(@{$self->templates}){
        for my $torsion(@{$template->fourmers}){
            push @lines, $torsion->string_repr;
        }
    }

    push @lines, "[Ramachandran data]\n";
    push @lines, $self->rama_data->string_repr;

    push @lines, "[Ramachandran]\n";
    push @lines, $self->query->ramachandran->string_repr;

    return join q{}, @lines;

}

sub TO_JSON {
    my ($self) = @_;

    my %json = (
        atoms   => [],
        linear  => [],
        angle   => [],
        torsion => [],
    );

    #Serialise the sequence into a flat string
    my @seq = map {$self->query->residues->{$_}}
        sort {$a <=> $b} keys %{$self->query->residues};
    $json{sequence} = join q{}, @seq;

    #Build the atom representation
    for my $res_i(sort {$a <=> $b} keys %{$self->query->backbone}){
        my $res = $self->query->backbone->{$res_i};
        push @{$json{atoms}}, $_ for @{$res->atoms};
    }

    #Build linear springs
    #Start with backbone springs
    push @{$json{linear}}, $_ for @{$self->query->backbone_springs};
    #Then add springs from templates
    for my $template(@{$self->templates}){
        my $pairs = $template->pairs;
        $pairs = $_->filter($pairs) for @{$self->spring_filters};
        push @{$json{linear}}, $_ for @{$pairs};
    }

    #Add angles from query
    push @{$json{angle}}, $_ for @{$self->query->angles};

    #Add torsion angles from each template
    for my $template(@{$self->templates}){
        push @{$json{torsion}}, $_ for @{$template->fourmers};
    }

    return \%json;
}


__PACKAGE__->meta->make_immutable;
1;
