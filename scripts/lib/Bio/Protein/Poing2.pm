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

=item C<until> (Default: C<synth_time * (length + 5)>

Run poing2 until the time reaches this value. The default is based on the
synthesis time and the length of the protein.

=cut

has until => (is => 'ro', isa => 'Num', lazy => 1, builder => '_build_until');

sub _build_until {
    my ($self) = @_;
    return $self->synth_time * ($self->query->length + 5);
}

=item C<synth_time> (Default: 500)

Time between each atom being synthesised.

=cut

has synth_time => (is => 'ro', isa => 'Num', default => 500);

=item C<drag_coefficient> (Default: -0.5)

Drag coefficient (when not using the shielded drag force).

=cut

has drag_coefficient => (is => 'ro', isa => 'Num', default => -0.5);

=item C<shield_drag> (Default: true)

Use the shielded drag force.

=cut

has shield_drag => (is => 'ro', isa => 'Bool', default => 1);

=item C<use_sterics> (Default: true)

Use the steric force

=cut

has use_sterics => (is => 'ro', isa => 'Bool', default => 1);

=item C<use_water> (Default: true)

Bombard the model with water molecules.

=cut

has use_water => (is => 'ro', isa => 'Bool', default => 1);

=item C<fix_before> (Default: 50, disabled if C<record_jitter> is false)

Requires this many atoms, from the most-recently synthesised atom backwards, to
remain free. The previous atoms may be fixed if their jitter becomes low enough.

=cut

has fix_before => (is => 'ro', isa => 'Int', default => 50);

=item C<record_time> (Default: 10 × timestep, disabled if C<record_jitter> is false)

Record the jitter at this time interval. Can be useful to avoid recording
millions of points when operating with a large C<fix_before> and and
C<synth_time> value.

=cut

has record_time => (is => 'ro', isa => 'Num', lazy => 1, builder => '_build_record_time');

sub _build_record_time {
    my ($self) = @_;
    return $self->timestep * 10;
}

=item C<max_jitter> (Default: 0.01 Å)

Atoms with average jitter below this value are frozen if the C<record_jitter>
option is set.

=cut

has max_jitter => (is => 'ro', isa => 'Num', default => 0.01);

=item C<record_jitter> (Default: false)

Should jitter be recorded and atoms near equilibrium be frozen? If true, the
C<fix_before> C<record_time> and C<max_jitter> options are stored in the
configuration.

=cut

has record_jitter => (is => 'ro', isa => 'Bool', default => 0);

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

=item C<verbose> [optional]

Print progress information to standard error.

=cut

has verbose => (is => 'ro', isa => 'Bool', default => 1);


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
        synth_time  => $self->synth_time + 0,
        until       => $self->until + 0,
        shield_drag => \$self->shield_drag,
        use_sterics => \$self->use_sterics,
        use_water   => \$self->use_water,
        drag_coefficient => $self->drag_coefficient + 0,
        atoms   => [],
        linear  => [],
        angle   => [],
        torsion => [],
        ramachandran => {},
    );
    if($self->record_jitter){
        $json{fix_before}  = $self->fix_before + 0;
        $json{record_time} = $self->record_time + 0;
        $json{max_jitter}  = $self->max_jitter + 0;
    }

    #Serialise the sequence into a flat string
    print STDERR "Building query sequence\n" if $self->verbose;
    my @seq = map {$self->query->residues->{$_}->oneletter}
        sort {$a <=> $b} keys %{$self->query->residues};
    $json{sequence} = join q{}, @seq;

    #Build the atom representation
    print STDERR "Building atom representation\n" if $self->verbose;
    for my $res_i(sort {$a <=> $b} keys %{$self->query->backbone}){
        my $res = $self->query->backbone->{$res_i};
        push @{$json{atoms}}, $_ for @{$res->atoms};
    }

    #Build linear springs
    #Start with backbone springs
    print STDERR "Building backbone springs\n" if $self->verbose;
    push @{$json{linear}}, $_ for @{$self->query->backbone_springs};
    #Then add springs from templates
    for my $template(@{$self->templates}){
        print STDERR "Building springs for ", $template->model, "\n"
            if $self->verbose;

        my $pairs = $template->pairs;
        $pairs = $_->filter($pairs) for @{$self->spring_filters};
        push @{$json{linear}}, $_ for @{$pairs};
    }

    #Add angles from query
    print STDERR "Building bond angles\n" if $self->verbose;
    push @{$json{angle}}, $_ for @{$self->query->angles};

    #Add torsion angles from each template
    for my $template(@{$self->templates}){
        print STDERR "Building torsions for ", $template->model, "\n"
            if $self->verbose;

        push @{$json{torsion}}, $_ for @{$template->fourmers};
    }

    #Ramachandran data and constraints
    print STDERR "Building Ramachandran constraints\n" if $self->verbose;
    $json{ramachandran}->{constraints} = $self->query->ramachandran;
    $json{ramachandran}->{data} = $self->rama_data;

    return \%json;
}


__PACKAGE__->meta->make_immutable;
1;
