package Bio::Protein::Poing2::Vector;
use strict;
use warnings;
use Carp;
use Bio::Protein::Poing2::Class;
use base 'Bio::Protein::Poing2::Class';

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

1;
