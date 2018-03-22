#!/usr/bin/env perl
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME => 'Bio::Protein::Poing2',
    VERSION_FROM => 'lib/Bio/Protein/Poing2.pm',
    TEST_REQUIRES  => {
        "Test::More" => 0,
        "Test::Exception" => 0,
    },
    PREREQ_PM => {
        "Try::Tiny" => 0,
        "Moose" => 0,
        "JSON::XS" => 0,
        "Carp" => 0,
        "Math::Vector::Real" => 0,
        "List::BinarySearch" => 0,
    },
    EXE_FILES => [qw(
        bin/build_config.pl
        bin/final_model.pl
        bin/pulchra_all.pl
    )],
);
