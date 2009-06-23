#!/usr/bin/perl

use strict;
use Module::Build;

my $SAMTOOLS = "$ENV{SAMTOOLS}";
die <<END unless -d $SAMTOOLS && -e "$SAMTOOLS/libbam.a";
Oops. Can\'t find libbam.a!

Please set the SAMTOOLS environment variable to point to a
samtools distribution directory containing a compiled libbam.a.
END
;

my $build = Module::Build->new(
    module_name        => 'Bio::SamTools',
    dist_version_from  => 'lib/Bio/DB/Sam.pm',
    dist_author        => 'Lincoln Stein <lincoln.stein@gmail.com>',
    dist_abstract      => 'The GMOD Generic Genome Browser',
    license            => 'perl',
    include_dirs       => [$SAMTOOLS],
    extra_linker_flags => ["-L$SAMTOOLS",'-lbam','-lz'],

    extra_compiler_flags=>[

        # must match DFLAGS in Samtools Makefile
	'-D_IOLIB=2','-D_FILE_OFFSET_BITS=64', 

	# turn off warnings originating in Perl's Newx* calls
	'-Wformat=0'
    ],

    build_requires => {
      'ExtUtils::CBuilder' => 0,
    },
    requires     => {
	'perl'                => '5.008',
    },
#    create_makefile_pl => 'passthrough',
    );

$build->create_build_script;
