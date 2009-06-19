package Bio::DB::Sam;
# $Id: Sam.pm,v 1.7 2009-06-19 20:21:38 lstein Exp $

=head1 NAME

Bio::DB::Sam -- Read SAM/BAM database files

=head1 SYNOPSIS

 use Bio::DB::Sam;

 # high level API
 my $sam = Bio::DB::Sam->new(-fasta=>"data/ex1.fa",
			     -bam  =>"data/ex1.bam");

 my @targets    = $sam->seq_ids;
 my @alignments = $sam->get_features_by_location(-seq_id => 'seq2',
                                                 -start  => 500,
                                                 -end    => 800);
 for my $a (@alignments) {
    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $cigar  = $a->cigar_str;
    my $paired = $a->get_tag_values('PAIRED');
    my $ref_dna   = $a->dna;        # reference sequence
    my $query_dna = $a->query->dna; # query sequence
    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match
 }

 my @pairs = $sam->get_features_by_location(-type   => 'read_pair',
                                            -seq_id => 'seq2',
                                            -start  => 500,
                                            -end    => 800);

 for my $pair (@pairs) {
    my $length                    = $pair->length;   # insert length
    my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
    my $f_start = $first_mate->start;
    my $s_start = $second_mate->start;
 }

 # low level API
 my $bam          = Bio::DB::Bam->open('/path/to/bamfile');
 my $header       = $bam->header;
 my $target_count = $header->n_targets;
 my $target_names = $header->target_name;
 while (my $align = $bam->read1) {
    my $seqid     = $target_names->[$align->tid];
    my $start     = $align->pos+1;
    my $end       = $align->calend;
    my $cigar     = $align->cigar_str;
 }

 my $index = Bio::DB::Bam->index_open('/path/to/bamfile');
 my $callback = sub {
     my $alignment = shift;
     my $start       = $alignment->start;
     my $end         = $alignment->end;
     my $seqid       = $target_names->[$alignment->tid];
     print $alignment->qname," aligns to $seqid:$start..$end\n";
 }
 my $header = $index->header;
 $index->fetch($bam,$header->parse_region('seq2'),$callback);

=head1 DESCRIPTION

This module provides a Perl interface to the libbam library for
indexed and unindexed SAM/BAM sequence alignment databases. It
provides support for retrieving information on individual alignments,
read pairs, and alignment coverage information across large
regions. It also provides callback functionality for calling SNPs and
performing other base-by-base functions. Most operations are
compatible with the BioPerl Bio::SeqFeatureI interface, allowing BAM
files to be used as a backend to the GBrowse genome browser
application (gmod.sourceforge.net).

=head2 The high-level API

The high-level API provides a BioPerl-compatible interface to indexed
BAM files. The BAM database is treated as a collection of
Bio::SeqFeatureI features, and can be searched for features by name,
location, type and combinations of feature tags such as whether the
alignment is part of a mate-pair.

When opening a BAM database using the high-level API, you provide the
pathnames of two files: the FASTA file that contains the reference
genome sequence, and the BAM file that contains the query sequences
and their alignments. If either of the two files needs to be indexed,
the indexing will happen automatically. You can then query the
database for alignment features by combinations of name, position,
type, and feature tag. The high-level API provides access to up to
four feature "types":

 * "match": The "raw" unpaired alignment between a read and the
   reference sequence.

 * "read_pair": Paired alignments; a single composite
   feature that contains two subfeatures for the alignments of each 
   of the mates in a mate pair.

 * "coverage": A feature that spans a region of interest that contains
   numeric information on the coverage of reads across the region.

 * "region": A way of retrieving information about the reference
   sequence. Searching for features of type "region" will return a
   list of chromosomes or contigs in the reference sequence, rather
   than read alignments.

 * "chromosome": A synonym for "region".

B<Features> can be en masse in a single call, retrieved in a
memory-efficient streaming basis using an iterator, or interrogated
using a filehandle that return a series of TAM-format lines.

B<SAM alignment flags> can be retrieved using BioPerl's feature "tag"
mechanism. For example, to interrogate the FIRST_MATE flag, one
fetches the "FIRST_MATE" tag:

  warn "aye aye captain!" if $alignment->get_tag_values('FIRST_MATE');

The Bio::SeqFeatureI interface has been extended to retrieve all flags
as a compact human-readable string, and to return the CIGAR alignment
in a variety of formats.  

B<Split alignments>, such as reads that cover introns, are dealt with
in one of two ways. The default is to leave split alignments alone:
they can be detected by one or more "N" operations in the CIGAR
string. Optionally, you can choose to have the API split these
alignments across two or more subfeatures; the CIGAR strings of these
split alignments will be adjusted accordingly.

B<Interface to the pileup routines> The API provides you with access
to the samtools "pileup" API. This gives you the ability to write a
callback that will be invoked on every column of the alignment for the
purpose of calculating coverage, quality score metrics, or SNP
calling.

The B<main object classes> that you will be dealing with in the
high-level API are as follows:

 * Bio::DB::Sam               -- A collection of alignments and reference sequences.
 * Bio::DB::Bam::Alignment    -- The alignment between a query and the reference.
 * Bio::DB::Bam::Query        -- An object corresponding to the query sequence.

You may encounter other classes as well. These include:

 * Bio::DB::Sam::Segment       -- This corresponds to a region on the reference
                                  sequence.
 * Bio::DB::Sam::Constants     -- This defines CIGAR symbol constants and flags.
 * Bio::DB::Bam::AlignWrapper  -- An alignment helper object that adds split
                                  alignment functionality. See Bio::DB::Bam::Alignment
                                  for the documentation on using it.
 * Bio::DB::Bam::ReadIterator  -- An iterator that mediates the one-feature-at-a-time 
                                  retrieval mechanism.
 * Bio::DB::Bam::FetchIterator -- Another iterator for feature-at-a-time retrieval.

=head2 The low-level API

The low-level API closely mirrors that of the libbam library. It
provides the ability to open TAM and BAM files, read and write to
them, build indexes, and perform searches across them. There is less
overhead to using the API because there is very little Perl memory
management, but the functions are less convenient to use. Some
operations, such as writing BAM files, are only available through the
low-level API.

The classes you will be interacting with in the low-level API are as
follows:

 * Bio::DB::Tam            -- Methods that read and write TAM (text SAM) files.
 * Bio::DB::Bam            -- Methods that read and write BAM (binary SAM) files.
 * Bio::DB::Bam::Header    -- Methods for manipulating the BAM file header.
 * Bio::DB::Bam::Index     -- Methods for retrieving data from indexed BAM files.
 * Bio::DB::Bam::Alignment -- Methods for manipulating alignment data.
 * Bio::DB::Bam::Pileup    -- Methods for manipulating the pileup data structure.
 * Bio::DB::Sam::Fai       -- Methods for creating and reading from indexed Fasta
                              files.

=head1 METHODS

We cover the high-level API first. The high-level API code can be
found in the files Bio/DB/Sam.pm, Bio/DB/Sam/*.pm, and
Bio/DB/Bam/*.pm.

=head2 Bio::DB::Sam Constructor and basic accessors

=over 4

=item $sam = Bio::DB::Sam->new(%options)

The Bio::DB::Sam object combines a Fasta file of the reference
sequences with a BAM file to allow for convenient retrieval of
human-readable sequence IDs and reference sequences. The new()
constructor accepts a -name=>value style list of options as
follows:

  Option         Description
  ------         -------------

  -fasta         Path to the Fasta file that contains
                   the reference sequences (required).

  -bam           Path to the BAM file that contains the
                   alignments (required).

  -expand_flags  A boolean value. If true then the standard
                   alignment flags will be broken out as 
                   individual tags such as 'UNMAPPED' (default
                   false).

  -split_splices A boolean value. If true, then alignments that
                  are split across splices will be broken out
                  into a single alignment containing two sub-
                  alignments (default false).

  -split          The same as -split_splices.

An example of a typical new() constructor invocation is:
 
  $sam = Bio::DB::Sam->new(-fasta => '/home/projects/genomes/hu17.fa',
                           -bam   => '/home/projects/alignments/ej88.bam',
                           -expand_flags  => 1,
                           -split_splices => 1);

B<-expand_flags> option, if true, has the effect of turning each of
the standard SAM flags into a separately retrievable B<tag> in the
Bio::SeqFeatureI interface. Otherwise, the standard flags will be
concatenated in easily parseable form as a tag named "FLAGS". See
get_all_tags() and get_tag_values() for more information.

Any two-letter extension flags, such as H0 or H1, will always appear
as separate tags regardless of the setting.

B<-split_splices> has the effect of breaking up alignments that
contain an "N" operation into subparts for more convenient
manipulation. For example, if you have both paired reads and spliced
alignments in the BAM file, the following code shows the subpart
relationships:

  $pair        = $sam->get_feature_by_name('E113:01:01:23');
  @mates       = $pair->get_SeqFeatures;
  @mate1_parts = $mates[0]->get_SeqFeatures;
  @mate2_parts = $mates[1]->get_SeqFeatures;

Because there is some overhead to splitting up the spliced alignments,
this option is false by default.

=item $flag = $sam->expand_flags([$new_value])

Get or set the expand_flags option. This can be done after object
creation and will have an immediate effect on all alignments fetched
from the BAM file.

=item $flag = $sam->split_splices([$new_value])

Get or set the split_splices option. This can be done after object
creation and will affect all alignments fetched from the BAM file
B<subsequently.>

=item $header = $sam->header

Return the Bio::DB::Bam::Header object associated with the BAM
file. You can manipulate the header using the low-level API.

=item $fai    = $sam->fai

Returns the Bio::DB::Sam::Fai object associated with the Fasta
file. You can then manipuate this object with the low-level API.

B<The index will be built automatically for you if it does not already
exist. If index building is necessarily, the process will need write
privileges to the same directory in which the Fasta file resides.> If
the process does not have write permission, then the call will fail.
Unfortunately, the BAM library does not do great error recovery for
this condition, and you may experience a core dump. This is not
trappable via an eval {}.

=item $bai    = $sam->bam_index

Return the Bio::DB::Bam::Index object associated with the BAM file. 

B<The BAM file index will be built automatically for you if it does
not already exist. In addition, if the BAM file is not already sorted
by chromosome and coordinate, it will be sorted automatically, an
operation that consumes significant time and disk space. The current
process must have write permission to the directory in which the BAM
file resides in order for this to work.> In case of a permissions
problem, the Perl library will catch the error and die. You can trap
it with an eval {}.

=back

=head2 Getting information about reference sequences

The Bio::DB::Sam object provides the following methods for getting
information about the reference sequence(s) contained in the
associated Fasta file.

=over 4

=item @seq_ids = $sam->seq_ids

Returns an unsorted list of the IDs of the reference sequences (known
elsewhere in this document as seq_ids). This is the same as the
identifier following the ">" sign in the Fasta file (e.g. "chr1").

=item $num_targets = $sam->n_targets

Return the number of reference sequences.

=item $length = $sam->length('seqid')

Returns the length of the reference sequence named "seqid".

=item $seq_id = $sam->target_name($tid)

Translates a numeric target ID (TID) returned by the low-level API
into a seq_id used by the high-level API.

=item $length = $sam->target_len($tid)

Translates a numeric target ID (TID) from the low-level API to a
sequence length.

=back

=head2 Creating and querying segments

The Bio::DB::Sam::Segment object refers to a region on the reference
sequence. It is used to retrieve the sequence of the reference, as
well as alignments that overlap with the region.

=over 4

=item $segment = $sam->segment($seqid,$start,$end);

=item $segment = $sam->segment(-seq_id=>'chr1',-start=>5000,-end=>6000);

Segments are created using the Bio:DB::Sam->segment() method. It can
be called using one to three positional arguments corresponding to the
seq_id of the reference sequence, and optionally the start and end
positions of a subregion on the sequence. If the start and/or end are
undefined, they will be replaced with the beginning and end of the
sequence respectively.

Alternatively, you may call segment() with named -seq_id, -start and
-end arguments.

All coordinates are 1-based.

=item $seqid = $segment->seq_id

Return the segment's sequence ID.

=item $start = $segment->start

Return the segment's start position.

=item $end  = $segment->end

Return the segment's end position.

=item $strand = $segment->strand

Return the strand of the segment (always 0).

=item $length = $segment->length

Return the length of the segment.

=item $dna    = $segment->dna

Return the DNA string for the reference sequence under this segment.

=item $seq    = $segment->seq

Return a Bio::PrimarySeq object corresponding to the sequence of the
reference under this segment. You can get the actual DNA string in
this redundant-looking way:

 $dna = $segment->seq->seq

The advantage of working with a Bio::PrimarySeq object is that you can
perform operations on it, including taking its reverse complement and
subsequences.

=item @alignments = $segment->features(%args)

Return alignments that overlap the segment in the associated BAM
file. The optional %args list allows you to filter features by name,
tag or other attributes. See the documentation of the
Bio::DB::Sam->features() method for the full list of options. Common
options include:

 # get all the overlapping alignments
 @all_alignments = $segment->features;  

 # get an iterator across the alignments
 my $iterator     = $segment->features(-iterator=>1);
 while (my $align = $iterator->next_seq) { do something }

 # get a TAM filehandle across the alignments
 my $fh           = $segment->features(-fh=>1);
 while (<$fh>) { print }

 # get only the alignments with unmapped mates
 my @unmapped    = $segment->features(-flags=>{M_UNMAPPED=>1});

 # get coverage across this region
 my ($coverage)       = $segment->features('coverage');
 my @data_points      = $coverage->coverage;

=item $tag = $segment->primary_tag

=item $tag = $segment->source_tag

Return the strings "region" and "sam/bam" respectively. These methods
allow the segment to be passed to BioPerl methods that expect
Bio::SeqFeatureI objects.

=item $segment->name, $segment->display_name, $segment->get_SeqFeatures, $segment->get_tag_values

These methods are provided for Bio::SeqFeatureI compatibility and
don't do anything of interest.

=back

=cut

use strict;
use warnings;

use Carp 'croak';
use Bio::SeqFeature::Lite;
use Bio::PrimarySeq;

use base 'DynaLoader';
our $VERSION = '0.03';
bootstrap Bio::DB::Sam;

use Bio::DB::Bam::Alignment;
use Bio::DB::Sam::Segment;
use Bio::DB::Bam::AlignWrapper;
use Bio::DB::Bam::FetchIterator;
use Bio::DB::Bam::ReadIterator;

sub new {
    my $class         = shift;
    my %args          = @_;
    my $fa_path       = $args{-fasta} or croak "-fasta argument required";
    my $bam_path      = $args{-bam}   or croak "-bam argument required";
    my $expand_flags  = $args{-expand_flags};
    my $split_splices = $args{-split} || $args{-split_splices};
    -e $fa_path  && -r _  or croak "$fa_path does not exist or is not readable";
    -e $bam_path && -r _  or croak "$fa_path does not exist or is not readable";
    my $fai = Bio::DB::Sam::Fai->open($fa_path)  or croak "$fa_path open: $!";
    my $bam = Bio::DB::Bam->open($bam_path)      or croak "$bam_path open: $!";
    my $self =  bless {
	fai           => $fai,
	bam           => $bam,
	bam_path      => $bam_path,
	expand_flags  => $expand_flags,
	split_splices => $split_splices,
    },ref $class || $class;
    $self->header;  # catch it
    return $self;
}

sub header {
    my $self = shift;
    return $self->{header} ||= $self->{bam}->header;
}

sub fai { shift->{fai} }

sub expand_flags {
    my $self = shift;
    my $d    = $self->{expand_flags};
    $self->{expand_flags} = shift if @_;
    $d;
}

sub split_splices {
    my $self = shift;
    my $d    = $self->{split_splices};
    $self->{split_splices} = shift if @_;
    $d;
}

sub reset_read {
    my $self = shift;
    $self->{bam}->header;
}

sub n_targets {
    shift->header->n_targets;
}

sub target_name {
    my $self = shift;
    my $tid  = shift;
    return $self->header->target_name->[$tid];
}

sub target_len {
    my $self = shift;
    my $tid  = shift;
    return $self->header->target_len->[$tid];
}

sub seq_ids {
    my $self    = shift;
    my $targets = $self->_cache_targets;
    return keys %{$targets};
}

sub _cache_targets {
    my $self = shift;
    return $self->{targets} if exists $self->{targets};
    my @targets = map {lc $_} @{$self->header->target_name};
    my @lengths =             @{$self->header->target_len};
    my %targets;
    @targets{@targets}      = @lengths;  # just you try to figure out what this is doing!
    return $self->{targets} = \%targets;
}


sub length {
    my $self        = shift;
    my $target_name = shift;
    return $self->_cache_targets->{lc $target_name};
}

sub _fetch {
    my $self     = shift;
    my $region   = shift;
    my $callback = shift;

    my $header              = $self->{bam}->header;
    my ($seqid,$start,$end) = $header->parse_region($region);
    return unless defined $seqid;
    my $index  = $self->bam_index;
    $index->fetch($self->{bam},$seqid,$start,$end,$callback,$self);
}

sub fetch {
    my $self     = shift;
    my $region   = shift;
    my $callback = shift;
    
    my $code     = sub {
	my ($align,$self) = @_;
	$callback->(Bio::DB::Bam::AlignWrapper->new($align,$self));
    };
    $self->_fetch($region,$code);
}

sub pileup {
    my $self   = shift;
    my ($region,$callback) = @_;

    my $header   = $self->header;
    my ($seqid,$start,$end) = $header->parse_region($region);
    return unless defined $seqid;

    my $refnames = $self->header->target_name;

    my $code = sub {
	my ($tid,$pos,$pileup) = @_;
	my $seqid = $refnames->[$tid];
	$callback->($seqid,$pos+1,$pileup);
    };

    my $index  = $self->bam_index;
    $index->pileup($self->{bam},$seqid,$start,$end,$code);
}

# segment returns a segment across the reference
# it will not work on a arbitrary aligned feature
sub segment {
    my $self = shift;
    my ($seqid,$start,$end) = @_;

    if ($_[0] =~ /^-/) {
	my %args = @_;
	$seqid = $args{-seq_id} || $args{-name};
	$start = $args{-start};
	$end   = $args{-stop}    || $args{-end};
    } else {
	($seqid,$start,$end) = @_;
    }

    my $targets = $self->_cache_targets;
    return unless exists $targets->{lc $seqid};

    $start = 1                     unless defined $start;
    $end   = $targets->{lc $seqid} unless defined $end;
    $start = 1 if $start < 1;
    $end   = $targets->{lc $seqid} if $end > $targets->{lc $seqid};

    return Bio::DB::Sam::Segment->new($self,$seqid,$start,$end);
}

sub get_features_by_location {
    my $self = shift;
    my %args;

    if ($_[0] =~ /^-/) { # named args
	%args = @_;
    } else {             # positional args
	$args{-seq_id} = shift;
	$args{-start}  = shift;
	$args{-end}    = shift;
    }
    $self->features(%args);
}

sub get_features_by_attribute {
  my $self       = shift;
  my %attributes = ref($_[0]) ? %{$_[0]} : @_;
  $self->features(-attributes=>\%attributes);
}

sub get_features_by_tag {
    shift->get_features_by_attribute(@_);
}

sub get_feature_by_name {
    my $self = shift;
    my %args;
    if ($_[0] =~ /^-/) {
	%args = @_;
    } else {
	$args{-name} = shift;
    }
    $self->features(%args);
}

sub get_features_by_name { shift->get_feature_by_name(@_) }

sub get_feature_by_id {
    my $self = shift;
    my $id   = shift;
    my ($name,$seqid,$start,$end,$strand) = map {s/%3B/;/ig;$_} split ';',$id;
    return unless $name && $seqid;
    my @features = $self->features(-name=>$name,
				   -seq_id=>$seqid,
				   -start=>$start,
				   -end=>$end,
				   -strand=>$strand);
    return unless @features;
    return $features[0];
}


sub get_seq_stream {
    my $self = shift;
    $self->features(-iterator=>1,@_);
}

sub get_seq_fh {
    my $self = shift;
    $self->features(-fh=>1,@_);
}

sub types {
    return qw(match read_pair coverage region);
}

sub features {
    my $self = shift;

    my %args;
    if ($_[0] !~ /^-/) {
	$args{-type} = \@_;
    } else {
	%args = @_;
    }

    my $seqid     = $args{-seq_id} || $args{-seqid};
    my $start     = $args{-start};
    my $end       = $args{-end}  || $args{-stop};
    my $types     = $args{-type} || $args{-types} || [];
    my $attributes = $args{-attributes} || $args{-tags} || $args{-flags};
    my $iterator  = $args{-iterator};
    my $fh        = $args{-fh};
    my $filter    = $args{-filter};

    $types        = [$types] unless ref $types;
    $types        = [$args{-class}] if !@$types && defined $args{-class};
    my $use_index = defined $seqid;

    # we do some special casing to retrieve target (reference) sequences
    # if they are requested
     if (defined($args{-name})
 	&& (!@$types || $types->[0]=~/region|chromosome/) 
	 && !defined $seqid) {
 	my @results = $self->_segment_search(lc $args{-name});
 	return @results if @results;
     } elsif (@$types && $types->[0] =~ /region|chromosome/) {
 	return map {$self->segment($_)} $self->seq_ids;
     }

    my %seenit;
    my @types = grep {!$seenit{$_}++} ref $types ? @$types : $types;
    @types    = 'match' unless @types;

    # the filter is intended to be inserted into a closure
    # it will return undef from the closure unless the filter
    # criteria are satisfied
    if (!$filter) {
	$filter = '';
	$filter   .= $self->_filter_by_name(lc $args{-name})
	    if defined $args{-name};
	$filter   .= $self->_filter_by_attribute($attributes)
	    if defined $attributes;
    }

    # Special cases for unmunged data
    if (@types == 1 && $types[0] =~ /^match/) {

	# if iterator is requested, and no indexing is possible,
	# then we directly iterate through the database using read1()
	if ($iterator && !$use_index) {
	    $self->reset_read;
	    my $code = eval "sub {my \$a=shift;$filter;1}";
	    die $@ if $@;
	    return Bio::DB::Bam::ReadIterator->new($self->{bam},$code);
	}

	# TAM filehandle retrieval is requested
	elsif ($fh) {
	    return $self->_features_fh($seqid,$start,$end,$filter);
	}

    }

    # otherwise we're going to do a little magic
    my ($features,@result);

    for my $t (@types) {

	if ($t =~ /^(match|read_pair)/) {
	    
	    # fetch the features if type is 'match' or 'read_pair'
	    $features = $self->_filter_features($seqid,$start,$end,$filter);

	    # for "match" just return the alignments
	    if ($t =~ /^(match)/) {
		push @result,@$features;
	    } 

	    # otherwise aggregate mate pairs into two-level features
	    elsif ($t =~ /^read_pair/) {
		$self->_build_mates($features,\@result);
	    }
	    next;
	}

	# create a coverage graph if type is 'coverage'
	# specify coverage:N, to create a map of N bins
	# units are coverage per bp
	# resulting array will be stored in the "coverage" attribute
	if ($t =~ /^coverage:?(\d*)/) {
	    my $bins = $1;
	    push @result,$self->_coverage($seqid,$start,$end,$bins,$filter);
	}
	
    }

    return $iterator ? Bio::DB::Bam::FetchIterator->new(\@result)
	             : @result;
}

sub _filter_features {
    my $self = shift;
    my ($seqid,$start,$end,$filter,$do_tam_fh) = @_;

    my @result;
    my $action = $do_tam_fh ? '\$self->header->view1($a)'
                            : 'push @result,Bio::DB::Bam::AlignWrapper->new($a,$self)';

    my $user_code;
    if (ref ($filter) eq 'CODE') {
	$user_code = $filter;
	$filter = '';
    }

    my $callback = defined($seqid) ? <<INDEXED : <<NONINDEXED;
sub {
    my \$a = shift;
    $filter
    return unless defined \$a->start;
    $action;
}
INDEXED
sub {
    my \$a    = shift;
    $filter
    $action;
}
NONINDEXED
    ;

    my $code = eval $callback;
    die $@ if $@;

    if ($user_code) {
	my $new_callback = sub {
	    my $a = shift;
	    $code->($a) if $user_code->($a);
	};
	$self->_features($seqid,$start,$end,$new_callback);
    } else {
	$self->_features($seqid,$start,$end,$code);
    }

    return \@result;
}

sub _features {
    my $self = shift;
    my ($seqid,$start,$end,$callback) = @_;

    if (defined $seqid) {
 	my $region = $seqid;
 	if (defined $start) { 
 	    $region   .= ":$start";
 	    $region   .= "-$end"   if defined $end;
 	}
 	$self->_fetch($region,$callback);
    } 

    else {
	$self->reset_read;
	while (my $b = $self->{bam}->read1) {
	    $callback->($b);
 	}
    }
}

# build mate pairs
sub _build_mates {
    my $self = shift;
    my ($src,$dest) = @_;

    my %read_pairs;
    for my $a (@$src) {
	my $name = $a->display_name;
	unless ($read_pairs{$name}) {
	    my $isize = $a->isize;
	    my $start = $isize >= 0 ? $a->start : $a->end+$isize+1;
	    my $end   = $isize <= 0 ? $a->end   : $a->start+$isize-1;
	    $read_pairs{$name} = 
		Bio::SeqFeature::Lite->new(
		    -display_name => $name,
		    -seq_id       => $a->seq_id,
		    -start => $start,
		    -end   => $end,
		    -type  => 'read_pair',
		    -class => 'read_pair',
		);
	}
	$read_pairs{$name}->add_SeqFeature($a);
    }
    push @$dest,values %read_pairs;
}

sub _coverage {
    my $self = shift;
    my ($seqid,$start,$end,$bins,$filter) = @_;

    # Currently filter is ignored. In reality, we should
    # turn filter into a callback and invoke it on each 
    # position in the pileup.
    croak "cannot calculate coverage unless a -seq_id is provided"
	unless defined $seqid;

    my $region = $seqid;
    if (defined $start) { 
	$region   .= ":$start";
	$region   .= "-$end"   if defined $end;
    }

    my $header     = $self->{bam}->header;
    my ($id,$s,$e) = $header->parse_region($region);
    return unless defined $id;

    # parse_region may return a very high value if no end specified
    $end   = $e >= 1<<29 ? $header->target_len->[$id] : $e;
    $start = $s+1;
    $bins ||= $end-$start+1;

    my $index      = $self->bam_index;
    my $coverage   = $index->coverage($self->{bam},
				      $id,$s,$e,
				      $bins);

    return Bio::SeqFeature::Coverage->new(
	-display_name => "$seqid coverage",
	-seq_id       => $seqid,
	-start        => $start,
	-end          => $end,
	-strand       => 0,
	-type         => "coverage:$bins",
	-class        => "coverage:$bins",
	-attributes   => { coverage => [$coverage] }
    );
}

sub _segment_search {
    my $self = shift;
    my $name = shift;

    my $targets = $self->_cache_targets;
    return $self->segment($name) if $targets->{$name};

    if (my $regexp = $self->_glob_match($name)) {
	my @results = grep {/^$regexp$/i} keys %$targets;
	return map {$self->segment($_)} @results;
    }

    return;
}

sub bam_index {
    my $self = shift;
    return Bio::DB::Bam->index($self->{bam_path});
}

sub _features_fh {
    my $self  = shift;
    my ($seqid,$start,$end,$filter) = @_;

    my $result = open my $fh,"-|";
    if (!$result) {  # in child
	$self->_filter_features($seqid,$start,$end,$filter,'do_fh'); # will print TAM to stdout
	exit 0;
    }
    return $fh;
    
}

sub tam_fh {
    my $self   = shift;
    return $self->features(-fh=>1);
}

# return a fragment of code that will be placed in the eval "" filter
# to eliminate alignments that don't match by name
sub _filter_by_name {
    my $self = shift;
    my $name = shift;

    my $frag = "my \$name=\$a->qname; defined \$name or return; ";

    if (my $regexp = $self->_glob_match($name)) {
	$frag .= "return unless \$name =~ /^$regexp\$/i;\n";
    } else {
	$frag .= "return unless lc \$name eq '$name';\n";
    }
}

# return a fragment of code that will be placed in the eval "" filter
# to eliminate alignments that don't match by attribute
sub _filter_by_attribute {
    my $self       = shift;
    my $attributes = shift;
    my $result;
    for my $tag (keys %$attributes) {
	$result .= "my \$value = lc \$a->get_tag_values('$tag');\n";
	$result .= "return unless defined \$value;\n";
	my @comps = ref $attributes->{$tag} eq 'ARRAY' 
	    ? @{$attributes->{$tag}} 
	    : $attributes->{$tag};
	my @matches;
	for my $c (@comps) {
	    if ($c =~ /^[+-]?[\deE.]+$/) { # numeric-looking argument
		push @matches,"\$value == $c";
	    }
	    elsif (my $regexp = $self->_glob_match($c)) {
		push @matches,"\$value =~ /^$regexp\$/i";
	    }
	    else {
		push @matches,"\$value eq lc '$c'";
	    }
	}
	$result .= "return unless " . join (' OR ',@matches) . ";\n";
    }
    return $result;
}

# turn a glob expression into a regexp
sub _glob_match {
    my $self = shift;
    my $term = shift;
    return unless $term =~ /(?:^|[^\\])[*?]/;
    $term =~ s/(^|[^\\])([+\[\]^{}\$|\(\).])/$1\\$2/g;
    $term =~ s/(^|[^\\])\*/$1.*/g;
    $term =~ s/(^|[^\\])\?/$1./g;
    return $term;
}

package Bio::SeqFeature::Coverage;

use base 'Bio::SeqFeature::Lite';

sub coverage {
    my $self       = shift;
    my ($coverage) = $self->get_tag_values('coverage');
    return wantarray ? @$coverage : $coverage;
}

package Bio::DB::Bam;

sub index {
    my $self = shift;
    my $path = shift;
    unless (-e "${path}.bai" && (-M $path >= -M "${path}.bai")) {
	# if bam file is not sorted, then index_build will exit.
	# we spawn a shell to intercept this eventuality
	print STDERR "[bam_index_build] creating index for $path\n" if -t STDOUT;

	my $result = open my $fh,"-|";
	die "Couldn't fork $!" unless defined $result;

	if ($result == 0) { # in child
	    # dup stderr to stdout so that we can intercept messages from library
	    open STDERR,">&STDOUT";  
	    $self->index_build($path);
	    exit 0;
	}

	my $mesg = <$fh>;
	$mesg  ||= '';
	close $fh;
	if ($mesg =~ /not sorted/i) {
	    print STDERR "[bam_index_build] sorting by coordinate...\n" if -t STDOUT;
	    $self->sort_core(0,$path,"$path.sorted");
	    rename "$path.sorted.bam",$path;
	    $self->index_build($path);
	} elsif ($mesg) {
	    die $mesg;
	}
    }
    return $self->index_open($path);
}


1;
__END__


=head1 EXAMPLES


=head1 SEE ALSO

L<Bio::Perl>

=head1 AUTHOR

Lincoln Stein E<lt>lincoln.stein@oicr.on.caE<gt>.
E<lt>lincoln.stein@bmail.comE<gt>

Copyright (c) 2009 Ontario Institute for Cancer Research.

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut

