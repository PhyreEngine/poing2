package Bio::Protein::Poing2::Modules;
use strict;
use warnings;

sub init_class_system {
    eval {
        require Moose;
        Moose->import;
    };
    return unless $@;

    eval {
        require Mouse;
        Mouse->import;
    };
    return unless $@;

    require parent;
    parent->import('Bio::Protein::Poing2::Class');
};

sub init_vector_lib {
    eval {
        require Math::Vector::Real;
        Math::Vector::Real->import
    };
    if($@){
        require Bio::Protein::Poing2::Vector;
        Bio::Protein::Poing2::Vector->import;
    }
}


#Load class syntax sugar
BEGIN {
    init_class_system();
    init_vector_lib();
};

1;
