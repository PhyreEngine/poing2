package Bio::Protein::Poing2::Vector;
use strict;
use warnings;
use Carp;

#Load class syntax sugar
BEGIN {
    if   (require Moose){ Moose->import }
    elsif(require Mouse){ Mouse->import }
    else {require parent; parent->import('Bio::Protein::Poing2::Class') }
};

has coords => (is => 'ro', default => sub{[0,0,0]});

use overload
    '-'   => 'minus',
    '+'   => 'add',
    '*'   => 'multiply',
    '/'   => 'divide',
    '=='  => 'equals',
    'neg' => 'unary_minus',
    '.'   => 'dot',
    'x'   => 'cross';

sub _verify_dimensions {
    my ($a, $b) = @_;
    my $dim_a = @{$a->coords};
    my $dim_b = @{$b->coords};
    croak "Dimensionality incorrect ($dim_a vs $dim_b)" if $dim_a != $dim_b;
}

sub equals {
    my ($self, $other) = @_;
    _verify_dimensions($self, $other);
    for(0 .. $#{$self->coords}){
        return 0 if $self->coords->[$_] != $other->coords->[$_];
    }
    return 1;
}

sub add {
    my ($self, $other) = @_;

    _verify_dimensions($self, $other);
    my $new_coords = [];
    for(0 .. $#{$self->coords}){
        push @{$new_coords}, $self->coords->[$_] + $other->coords->[$_];
    }
    return Bio::Protein::Poing2::Vector->new(coords => $new_coords);
}

sub multiply {
    my ($self, $other) = @_;
    if(!ref($other)){
        return $self->_multiply_scalar($other);
    }else{
        croak "Vector multiplication not implemented.";
    }
}

sub _multiply_scalar {
    my ($self, $scalar) = @_;
    my $new_coords = [];
    for(0 .. $#{$self->coords}){
        push @{$new_coords}, $self->coords->[$_] * $scalar;
    }
    return Bio::Protein::Poing2::Vector->new(coords => $new_coords);
}

sub divide {
    my ($self, $other) = @_;
    return $self->multiply(1 / $other);
}

sub minus {
    my ($self, $other, $swap) = @_;
    my $new = $self->add(-$other);
    $new = -$new if $swap;
    return $new;
}

sub unary_minus {
    my ($self) = @_;
    return $self * -1;
}

sub dot {
    my ($self, $other) = @_;
    _verify_dimensions($self, $other);

    my $sum = 0;
    for(0 .. $#{$self->coords}){
        $sum += $self->coords->[$_] * $other->coords->[$_];
    }
    return $sum;
}

sub cross {
    my ($self, $other) = @_;
    croak "Vector product only makes sense in 3D" unless @{$self->coords} == 3;

    my $a1 = $self->coords;
    my $a2 = $other->coords;
    my $cross = [
        $a1->[1]*$a2->[2] - $a1->[2]*$a2->[1],
        $a1->[2]*$a2->[0] - $a1->[0]*$a2->[2],
        $a1->[0]*$a2->[1] - $a1->[1]*$a2->[0],
    ];
    return Bio::Protein::Poing2::Vector->new(coords => $cross);
}

sub mag {
    my ($self) = @_;
    return sqrt($self . $self);
}

#Rotate self around a given axis
sub vrot_axis {
    my ($self, $axis, $theta) = @_;

    my $ux = $axis->coords->[0];
    my $uy = $axis->coords->[1];
    my $uz = $axis->coords->[2];

    my $vx = $self->coords->[0];
    my $vy = $self->coords->[1];
    my $vz = $self->coords->[2];

    my $ct = cos($theta);
    my $st = sin($theta);

    my $coords = [0, 0, 0];
    $coords->[0] = $vx * ($ct + $ux*$ux*(1-$ct))
        + $vy * ($ux*$uy*(1-$ct) - $uz*$st)
        + $vz * ($ux*$uz*(1-$ct) + $uy*$st);

    $coords->[1] = $vx * ($uy*$ux*(1-$ct) + $uz*$st)
        + $vy * ($ct + $uy*$uy*(1-$ct))
        + $vz * ($uy*$uz*(1-$ct) - $ux*$st);

    $coords->[2] = $vx * ($uz*$ux*(1-$ct) - $uy*$st)
        + $vy * ($uz*$uy*(1-$ct) + $ux*$st)
        + $vz * ($ct + $uz*$uz*(1-$ct));

    return ref($self)->new(coords => $coords);
}

if(defined __PACKAGE__->meta){
    __PACKAGE__->meta->make_immutable;
}

1;
