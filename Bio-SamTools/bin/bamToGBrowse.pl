#!/usr/bin/perl

use strict;
use Bio::DB::Sam;

my $has_bigwig = eval {require Bio::DB::BigFile;1};

@ARGV >= 1 or die <<USAGE;
Usage: $0 <directory_path> [<ref.fa>]

Given the path to a directory and the fasta file for the reference
sequence, do the following:

 1. Find all SAM files in the indicated directory and convert them
    into BAM. These must end in one of the extensions ".sam" or ".sam.gz".
    A series of <base>.bam files will be created.

 2. Sort the newly created BAM files.

 3. Index BAM files that need indexing. This step will look for
      files named <base>_sorted.bam

 4. Create a set of BigWig files representing the coverage graph. These
      will be named <base>.bw.

 5. Create a skeletal GBrowse config file named "gbrowse.conf" that
    serves as a starting point for viewing these files. Previous versions
    of this file will be appended to.

If the Fasta file is not provided, then this script will look in the
designated directory for ONE .fa file to use.
USAGE
    ;

my($dir,$fasta) = @ARGV;
my $converter = BamToGBrowse->new($dir,$fasta,'verbose');
$converter->run();

exit 0;

package BamToGBrowse;
use Carp 'croak';
use File::Spec;
use File::Basename 'basename';
use File::Temp;

sub new {
    my $class = shift;
    my ($dir,$fasta,$debug) = @_;
    $dir or croak "Usage: $class->new(\$dir_to_index,[\$fasta_path])";
    $fasta ||= $class->find_fasta($dir) or croak "Cannot find a suitable fasta file in $dir";
    return bless { dir   => $dir,
		   fasta => $fasta,
		   debug => $debug||0
    },ref $class || $class;
}

sub run {
    my $self = shift;
    $self->sam_to_bam;
    $self->index_bam;
    $self->bam_to_wig;
    $self->make_conf;
}

sub dir   {shift->{dir}   }
sub fasta {shift->{fasta} }
sub debug {shift->{debug} }
sub msg {
    my $self = shift;
    return unless $self->debug;
    print STDERR @_,"\n";
}
sub err {
    my $self = shift;
    print STDERR @_,"\n";
}
sub files {
    my $self = shift;
    my @extensions = @_; # e.g. '.sam','.sam.gz';
    my $dir = $self->dir;
    return map {glob($self->dir_path("*$_"))} @extensions;
}
sub mtime {
    my $self = shift;
    my $file = shift;
    my $mtime = (stat($file))[9];
    return $mtime;
}
sub up_to_date {
    my $self = shift;
    my ($source,$target) = @_;
    return unless -e $target;
    return unless $self->mtime($target) > $self->mtime($source);
    return 1;
}
sub find_fasta {
    my $self  = shift;
    my $dir   = shift;
    my @files = glob(File::Spec->catfile($dir,"*.fa"));
    return unless @files == 1;
    return $files[0];
}
sub sam_to_bam {
    my $self = shift;
    $self->msg('Searching for SAM files');
    my @sam = $self->files($self->sam_extensions);

    $self->msg('Found ',@sam+0,' sam files');
    $self->convert_one_sam($_) foreach @sam;
}

sub sam_extensions { return qw(.sam .sam.gz) }

sub index_bam {
    my $self = shift;
    $self->msg('Searching for BAM files');
    my @bam = $self->files('.bam');
    $self->msg('Found ',@bam+0,' bam files');
    $self->index_one_bam($_) foreach @bam;
}

sub convert_one_sam {
    my $self = shift;
    my $sam  = shift;
    my $base = basename($sam,$self->sam_extensions);
    my $bam    = $self->dir_path("$base.bam");
    my $sorted = $self->dir_path("${base}_sorted.bam");

    if ($self->up_to_date($sam,$bam) or $self->up_to_date($sam,$sorted)) {
	$self->msg("$bam: up to date");
	return;
    }

    $self->msg("$bam: creating");

    # This is to create the .fai file. Do this in a block so that handle
    # goes out of scope when not needed.
    my $fasta = $self->fasta;
    -r $fasta or croak "$fasta is not readable";
    my $fai  = $fasta.".fai";

    unless ($self->up_to_date($fasta,$fai))
    {
	my $fai = Bio::DB::Sam::Fai->load($fasta)
	    or die "Could not load reference FASTA file for indexing this SAM file: $!";
    }

    my $tam = Bio::DB::Tam->open($sam)
	or die "Could not open SAM file for reading: $!";

    my $header = $tam->header_read2($fai);

    my $out = Bio::DB::Bam->open($bam,'w')
	or die "Could not open BAM file for writing: $!";

    $out->header_write($header);
    my $alignment = Bio::DB::Bam::Alignment->new();
    my $lines = 0;

    while ($tam->read1($header,$alignment) > 0) {
	$out->write1($alignment);
	$self->msg("converted $lines lines...") if ++$lines%100000 == 0;
    }
    undef $tam;
    undef $out; 

    $self->msg("converted $lines lines");
    $self->sort_bam($bam);
}

sub sort_bam {
    my $self = shift;
    my $bam  = shift;
    $self->msg("sorting $bam");
    my $basename = basename($bam,'.bam');

    my $sorted = $self->dir_path($basename.'_sorted');
    Bio::DB::Bam->sort_core(0,$bam,$sorted);
    unlink $bam;  # we don't want the unsorted version!
    return $sorted.'.bam';
}

# This guy is a little tricky because unsorted BAM files
# will terminate the process. We try to get around this by
# forking and reading STDERR to see what happened (the Samtools
# library is not great at returning status codes.
sub index_one_bam {
    my $self = shift;
    my $bam  = shift;
    my $base = basename($bam,'.bam');
    my $idx  = "${base}.bam.bai";

    if ($self->up_to_date($bam,$idx)) {
	$self->msg("$bam: index is up to date");
	return;
    }

    $self->msg("indexing $bam");
    my $err = $self->_fork_and_index_bam($bam);

    if ($err =~ /alignment is not sorted/) {
	$self->msg("$bam needs sorting");
	$bam = $self->sort_bam($bam);
	$err = $self->_fork_and_index_bam($bam);
    }
    
    if ($err) {
	$self->err("Could not index $bam: $err");
    } else {
	$self->msg("$bam indexed successfully");
    }
}

sub _fork_and_index_bam {
    my $self = shift;
    my $bam  = shift;
    my $pid = open my $pipe,"-|";
    die "Couldn't fork: $!" unless defined $pid;

    if ($pid) { # I am the parent
	my $output = join '',<$pipe>;
	return $output;
    }

    # Otherwise, I am the child
    open STDERR,">&STDOUT"; # get stderr to go to the pipe
    Bio::DB::Bam->index_build($bam);
    exit 0;
}
sub dir_path {
    my $self = shift;
    my $filename = shift;
    return File::Spec->catfile($self->dir,$filename);
}
sub bam_to_wig {
    my $self  = shift;
    $self->msg('Searching for .bai files');
    my @files = map {$self->dir_path(basename($_,'.bai'))} $self->files('.bai');
    $self->msg('Found ', @files+0,' files');
    $self->wiggle_one_bam($_) foreach @files;
}

sub wiggle_one_bam {
    my $self = shift;
    my $bam  = shift;

    my $base        = basename($bam,'.bam');
    my $bigwig      = $self->dir_path($base.'.bw');
    if ($self->up_to_date($bam,$bigwig)) {
	$self->msg("$bigwig is up to date");
	return;
    }

    if (-r '/dev/stdin' && -c _) {  # only works with linux, I think
	$self->_wiggle_one_bam_pipe($bam,$bigwig);
    } else {
	$self->_wiggle_one_bam_tempfile($bam,$bigwig);
    }
}

sub _wiggle_one_bam_pipe {
    my $self = shift;
    my ($bam,$bigwig)  = @_;

    my $pid  = open my $pipe,"|-";
    defined $pid or die "Couldn't fork: $!";

    if ($pid) { # I'm the parent; my job is to write the coverage data to stdout
	$self->write_coverage($bam,$pipe);
	close $pipe;
	return;
    }

    else {   # I'm the child; my job is to create the BigWig file from /dev/stdin
	$self->msg("writing bigwig file");
	my $chrom_sizes = $self->fasta.".fai";
	Bio::DB::BigFile->createBigWig('/dev/stdin',$chrom_sizes,$bigwig);
	exit 0;
    }
}

sub write_coverage {
    my $self = shift;
    my ($bamfile,$fh) = @_;

    $self->msg("calculating coverage for $bamfile");
    my $bam = Bio::DB::Sam->new(-bam=>$bamfile) 
	or die "Couldn't open $bamfile: $!";

    my $callback = sub {
	my ($seqid,$pos,$pileup,$sam) = @_;
	print $fh $pos,"\t",scalar @$pileup,"\n";
    };

    for my $seq_id ($bam->seq_ids) {
	print $fh "variableStep chrom=$seq_id span=1\n";
	$bam->fast_pileup($seq_id,$callback);
    }
}

sub _wiggle_one_bam_tempfile {
    my $self = shift;
    my ($bam,$bigwig) = @_;
    my $tmpfh = File::Temp->new(TEMPLATE => 'wigfileXXXXX',
				UNLINK   => 1,
				DIR      => $self->dir,
				SUFFIX   => '.wig');
    $self->write_coverage($bam,$tmpfh);
    close $tmpfh;

    $self->msg("writing bigwig file");
    my $chrom_sizes = $self->fasta.".fai";
    Bio::DB::BigFile->createBigWig($tmpfh,$chrom_sizes,$bigwig);
}

sub make_conf {
    my $self = shift;
    my $conf = $self->dir_path('gbrowse.conf');

    my $existing_config = -e $conf ? $self->parse_conf($conf) : {};

    my %tracks  = map  { basename($_,'.bw')=>1   } $self->files('.bw');
    my @new     = grep { !$existing_config->{$_} } keys %tracks;

    my $newfile = "$conf.new";
    open my $f,">",$newfile or die "$newfile: $!";

    for my $track (sort keys %tracks) {
	if ($existing_config->{$track}) {
	    print $f $existing_config->{$track};
	} else {
	    print $f $self->make_gbrowse_conf($track);
	}
    }
    close $f;
    rename $newfile,$conf;
}

sub parse_conf {
    my $self = shift;
    my $conf = shift;

    open my $f,$conf or die "$conf: $!";
    my ($current,%data);
    while (<$f>) {
	if (/^\[([^:]+)/) {
	    $current = $1;
	    $data{$current} = $_;
	} elsif ($current) {
	    $data{$current} .= $_;
	}
    }
    return \%data;
}

sub make_gbrowse_conf {
    my $self  = shift;
    my $track = shift;
    $self->msg("creating gbrowse stanza for $track");

    my $fasta = File::Spec->rel2abs($self->fasta);
    my $bam   = File::Spec->rel2abs($self->dir_path("$track.bam"));
    my $bw    = File::Spec->rel2abs($self->dir_path("$track.bw"));
    (my $key = $track) =~ s/_sorted//;

    my $result = <<END;
[$track:database]
db_adaptor = Bio::DB::Sam
db_args    = -fasta $fasta
      	     -bam   $bam
	     -split_splices 1

[${track}_bw:database]
db_adaptor = Bio::DB::BigWig
db_args    = sub { require Bio::DB::Sam;
                   return ( 
                       -bigwig => '$bw',
		       -fasta  => Bio::DB::Sam::Fai->open('$fasta'),
		       );
                 }
                        

[$track]
database = $track
feature  = read_pair
glyph    = segments
draw_target  = 1
show_mismatch= 1
mismatch_only = 1
mismatch_color = orange
indel_color    = yellow
bgcolor      = black
fgcolor      = black
height       = 4
label        = 1
label_position = left
label density = 50
bump         = fast
connector    = sub {
		  my \$glyph = pop;
		  return \$glyph->level == 0 ? 'dashed' : 'solid';
               }
maxdepth     = 2
box_subparts = 2
key          = Reads from $key

[$track:2000]
database = ${track}_bw
feature  = summary
glyph    = wiggle_whiskers

END
}
