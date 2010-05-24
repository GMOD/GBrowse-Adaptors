#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Bio::DB::Sam::SamToGBrowse;

my $has_bigwig = eval {require Bio::DB::BigFile;1};

@ARGV >= 1 or die <<USAGE;
Usage: $0 <directory_path> [<ref.fa>]

Given the path to a directory and the fasta file for the reference
sequence, do the following:

 1. Find all SAM files in the indicated directory and convert them
    into BAM. These must end in one of the extensions ".sam" or ".sam.gz".
    A series of <base>.bam files will be created.

 2. Sort the newly created BAM files, creating a series of files named
    <base>_sorted.bam.

 3. Index BAM files that need indexing. This step will look for
      files named <base>_sorted.bam

 4. Create a set of BigWig files representing the coverage graph. These
      will be named <base>.bw.

 5. Create a skeletal GBrowse config file named "gbrowse.conf" that
    serves as a starting point for viewing these files. Previous versions
    of this file will be appended to.

You can prepopulate the directory with sorted and indexed BAM files,
in which case the script will not attempt to resort or reindex them.
Unless the Fasta file is explicitly provided, this script will look in
the designated directory for ONE .fa file to use.

Note that you will need temporary space in the directory equivalent to
the size of the largest SAM file being processed. In addition, you
should have a copy of the bedGraphToBigWig executable somewhere on
your path (obtainable in binary form from
http://hgdownload.cse.ucsc.edu/admin/exe or source form from
http;//hgdownload.cse.ucsc.edu/admin/jksrc.zip). If the binary is not
found, this script will use its built-in BigWig loader, which uses a
lot of RAM.
USAGE
    ;

my($dir,$fasta) = @ARGV;
my $converter = Bio::DB::Sam::SamToGBrowse->new($dir,$fasta,'verbose');
$converter->run();

exit 0;

