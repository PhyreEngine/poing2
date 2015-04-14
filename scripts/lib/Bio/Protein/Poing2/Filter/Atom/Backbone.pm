package Bio::Protein::Poing2::Filter::Atom::Backbone;
use strict;
use warnings;
use Bio::Protein::Poing2::Data;

#Load class syntax sugar
BEGIN {
    if   (require Moose){ Moose->import }
    elsif(require Mouse){ Mouse->import }
    else {require parent; parent->import('Bio::Protein::Poing2::Class') }
};

=head1 NAME

Bio::Protein::Poing2::Filter::Atom::Backbone - Keep backbone atoms only

=head1 SYNOPSIS

    my $filt = Bio::Protein::Poing2::Filter::Atom::Backbone->new();
    my $filtered_residues = $filt->filter($hashref_of_residues);

=cut

sub filter {
    my ($self, $residues) = @_;
    for my $res_num(keys %{$residues}){
        my $res = $residues->{$res_num};
        $res->atoms([grep {
            $Bio::Protein::Poing2::Data::backbone{$_->name}
        } @{$res->atoms}]);
    }
    return $residues;
}

if(defined __PACKAGE__->meta){
    __PACKAGE__->meta->make_immutable;
}
1;
