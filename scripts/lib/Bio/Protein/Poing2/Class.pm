package Bio::Protein::Poing2::Class;
use strict;
use warnings;
use utf8;
use Carp;
use Exporter;

#We're just following Moose here, so ignore Moose
##no critic (Modules::ProhibitAutomaticExportation)
use base 'Exporter';
our @EXPORT = qw(has);
##use critic

=head1 NAME

Bio::Protein::Poing2::Class - Provides boilerplate for OO code

=head1 SYNOPSIS

    #In your package
    package Bio::Protein::Poing2::Foo;
    use parent 'Bio::Protein::Poing2::Class';

    has name   => (is => 'ro', require => 1);
    has coords => (is => 'rw', default => sub { [0, 0, 0] });

    #In a caller package

    #Dies: name is mandatory
    Bio::Protein::Poing2::Foo->new();

    my $obj = Bio::Protein::Poing2::Foo->new(name => 'C');
    print $obj->coords->[0]; #0
    $obj->coords([1, 2, 3]);

=head1 DESCRIPTION

This module is a very basic replacement for L<Moose>. I don't like to reinvent
the wheel, but there is often some resistance to installing dependency-heavy
modules like Moose on academic machines. This provides a small replacement that
takes the hassle away from defining constructors and attributes.

Note that this is probably a bad idea. It's too clever by half, and reinvents a
(probably slightly broken) tiny subset of Moose. With that said, it'll serve my
purposes.

=head1 METHODS

=head2 C<has($attr, %list)>

This method behaves similarly to Moose method with the same name. That is, it
creates getters (and optionally setters) for the attribute C<$attr>. The
attribute definition is stored in the package variable C<$meta> of the calling
package.

=cut

#Just stuff our attribute into $pkg->{meta}{attributes} so AUTOLOAD can find it

#Briefly disable Perl::Critic; we're just following Moose
##no critic (Subroutines::ProhibitSubroutinePrototypes)
sub has(@) {
##use critic
    my ($attr, %list) = @_;
    croak "Only call `has' in void context." if defined wantarray;

    my ($pkg) = caller;
    my $meta = do {
        ##no critic (TestingAndDebugging::ProhibitNoStrict)
        no strict 'refs';
        ${"${pkg}::meta"} ||= {};
        ${"${pkg}::meta"}
        ##use critic
    };
    $meta->{attributes} ||= {};

    croak "Attribute `$attr' already defined"
            if exists $meta->{attributes}{$attr};

    $meta->{attributes}{$attr} = \%list;

    if($list{is} eq 'rw'){
        ##no critic (TestingAndDebugging::ProhibitNoStrict)
        no strict 'refs';
        *{"${pkg}::$attr"} = sub {
            my ($self, $value) = @_;
            $self->{$attr} = $value if @_ == 2;
            return $self->{$attr};
        };
        ##use critic
    }else{
        ##no critic (TestingAndDebugging::ProhibitNoStrict)
        no strict 'refs';
        *{"${pkg}::$attr"} = sub {
            my ($self, $value) = @_;
            croak "Attribute `$attr' is read-only." if @_ > 1;
            return $self->{$attr};
        };
        ##use critic
    }
    return;
}


=head2 C<new(%list)>

Inherited by calling packages, this method implements the required argument
parsing for packages.

=cut

sub new {
    my ($pkg, %list) = @_;

    my %to_bless = ();

    my $meta = do {
        ##no critic (TestingAndDebugging::ProhibitNoStrict)
        no strict 'refs';
        ${"${pkg}::meta"}
        ##use critic
    };


    #Collect attributes from meta key. If no atts have been set, just bless the
    #hashref and return it.
    if($meta && $meta->{attributes}){

        #Go through all defined attributes and place default values in the
        #hash, or die on mandatory args that have not been supplied.
        my $atts = $meta->{attributes};
        for my $att(keys %{$atts}){

            #Check whether a required value is supplied
            if($atts->{$att}{required} && !$list{$att}){
                croak "Mandatory argument `$att' not supplied";
            }

            #Place the attribute into the hash with the default value. This can
            #be a code ref, so eval it if it is
            my $default = $atts->{$att}{default};
            if($default && ref($default) eq 'CODE'){
                $default = eval {&$default};
                if(!$default){
                    croak "Error evaluating default value for attribute `$att': $@";
                }
            }
            $to_bless{$att} = $default || undef;
        }

        #Go through the supplied arguments and make sure they are allowed,
        #setting values if supplied.
        for my $att(keys %list){
            croak "Unknown attribute `$att'" if !$atts->{$att};
            $to_bless{$att} = $list{$att};
        }
    }else{
        %to_bless = %list;
    }

    return bless \%to_bless, $pkg;
}

1;

