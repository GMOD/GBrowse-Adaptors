package Bio::DB::BigBed;

#$Id$

use strict;
use warnings;
use base 'Bio::DB::BigWig';
use Carp 'croak';

# high level interface to BigBed files

sub new {
    my $self = shift;
    my %args = @_;

    my $bw_path       = $args{-bigbed}  or croak "-bigbed argument required";
    my $fa_path       = $args{-fasta};
    my $dna_accessor  = $self->new_dna_accessor($fa_path);
    
    unless ($self->is_remote($bw_path)) {
	-e $bw_path or croak "$bw_path does not exist";
	-r _  or croak "is not readable";
    }

    my $bw = Bio::DB::BigFile->bigBedFileOpen($bw_path)
	or croak "$bw_path open: $!";

    return bless {
	bw => $bw,
	fa => $dna_accessor
    },ref $self || $self;
}

sub bb     { shift->bw  }
sub bigbed { shift->bw }

sub _type_to_iterator {
    my $self = shift;
    my $type = shift;

    return 'Bio::DB::BigBed::BinIterator'      if $type =~ /^bin/;
    return 'Bio::DB::BigBed::SummaryIterator'  if $type =~ /^summary/;
    return 'Bio::DB::BigBed::FeatureIterator'  if $type =~ /^(region|feature)/;
    return 'Bio::DB::BigBed::EmptyIterator';
}

############################################################

package Bio::DB::BigBed::FeatureIterator;

use base 'Bio::DB::BigWig::IntervalIterator';

sub _query {
    return 'bigBedInterval';
}

sub _make_feature {
}

############################################################

package Bio::DB::BigWig::BinIterator;

use Carp 'croak';

sub new {
    my $self   = shift;
    my ($bigwig,$options) = @_;
    my $bw     = $bigwig->bw;

    my (undef,$bins) = $options->{-type} =~ /^(bin)(?::(\d+))?/i
	or croak "invalid call to _get_bin_stream. -type argument must be bin[:bins]";
    $bins ||= 1024;

    my $arry   = $bw->bigWigSummaryArrayExtended($options->{-seq_id},
						 $options->{-start}-1,
						 $options->{-end},
						 $bins)
	or return;

    my $chrom_end = $bw->chromSize($options->{-seq_id});
    my $binsize   = ($options->{-end}-$options->{-start}+1)/$bins;
    return bless {
	array   => $arry,
	bigwig  => $bigwig,
	start   => $options->{-start},
	end     => $chrom_end,
	binsize => $binsize,
	seq_id  => $options->{-seq_id},
    },ref $self || $self;
}

sub next_seq {
    my $self = shift;
    my $filter = $self->{options}{-filter};

    my $array  = $self->{array};
    my $i      = shift @$array;
    my $f;
    while ($i) {
	my $end = int($self->{start}+$self->{binsize});
	$end    = $self->{end} if $end > $self->{end};
	$f = Bio::DB::BigWig::Feature->new(-seq_id => $self->{seq_id},
					   -start  => int($self->{start}),
					   -end    => $end,
					   -score  => $i,
					   -type   => 'bin',
					   -fa     => $self->{bigwig}->fa,
	    );
	$self->{start} += $self->{binsize} + 1;
	last if !$filter || $filter->($f);
	$i = shift @$array;
    }

    return $f if $i;
    return;
}

##################################################################
package Bio::DB::BigWig::SummaryIterator;

sub new {
    my $self = shift;
    my ($bigwig,$options) = @_;
    my $bw     = $bigwig->bw;
    return bless {
	feature => Bio::DB::BigWig::Summary->new(-seq_id => $options->{-seq_id},
						 -start  => $options->{-start},
						 -end    => $options->{-end},
						 -type   => 'summary',
						 -fa     => $bigwig->fa,
						 -bw     => $bigwig)
    },ref $self || $self;
}

sub next_seq {
    my $self = shift;
    my $d = $self->{feature};
    $self->{feature} = undef;
    return $d;
}

##################################################################

package Bio::DB::BigWig::EmptyIterator;

sub new { my $self = shift; return bless {},ref $self || $self }
sub next_seq { return }


##################################################################

package Bio::DB::BigWig::Feature;

use base 'Bio::Graphics::Feature';

sub new {
    my $self = shift;
    my $feat = $self->SUPER::new(@_);
    my %args = @_;
    $feat->{fa} = $args{-fa} if $args{-fa};
    return $feat;
}

sub dna {
    my $self = shift;
    my $fa     = $self->{fa} or return 'N' x $self->length;
    my $seq_id = $self->seq_id;
    my $start  = $self->start;
    my $end    = $self->end;
    return $fa->seq($seq_id,$start,$end);
}

sub seq {
    my $self = shift;
    return Bio::PrimarySeq->new(-seq=>$self->dna);
}


##################################################################

package Bio::DB::BigWig::Summary;

use base 'Bio::DB::BigWig::Feature';

sub new {
    my $self = shift;
    my $feat = $self->SUPER::new(@_);
    my %args = @_;
    $feat->{bw} = $args{-bw} if $args{-bw};
    return $feat;
}

sub statistical_summary {
    my $self = shift;
    my $bins = shift;
    $bins ||= 1024;

    my $bw = $self->{bw}->bw or return;
    return $bw->bigWigSummaryArrayExtended($self->seq_id,
					   $self->start-1,
					   $self->end,
					   $bins);
}

##################################################################
package Bio::DB::BigWig::Segment;

use base 'Bio::DB::BigWig::Summary';

sub features {
    my $self = shift;
    return $self->{bw}->features(-seq_id => $self->seq_id,
				 -start  => $self->start,
				 -end    => $self->end,
				 -type   => $_[0] || 'region');
}

sub get_seq_stream {
    my $self = shift;
    return $self->{bw}->get_seq_stream(-seq_id => $self->seq_id,
				       -start  => $self->start,
				       -end    => $self->end,
				       -type   => $_[0] || 'region');
}

##################################################################

package Bio::DB::BigWig::GlobalIterator;

sub new {
    my $self = shift;
    my ($bigwig,$inner_iterator,$options) = @_;
    my $cl = $bigwig->bw->chromList or return;

    my $s =  bless {
	cl_head  => $cl,     # keep in scope so not garbage collected
	current  => $cl->head,
	bigwig   => $bigwig,
	options  => $options,
	inner_i  => $inner_iterator,
    },ref $self || $self;

    $s->{interval}  = $s->_new_interval($bigwig,$cl->head,$options);
    return $s;
}


sub next_seq {
    my $self = shift;
    my $c    = $self->{current} or return;

    my $next = $self->{interval}->next_seq;
    return $next if $next;
    
    # if we get here, then there are no more intervals on current chromosome
    # try more chromosomes
    while (1) {
	$self->{current}  = $self->{current}->next or return;  # out of chromosomes
	$self->{interval} = $self->_new_interval($self->{bigwig},$self->{current},$self->{options});
	my $next = $self->{interval}->next_seq;
	return $next if $next;
    }
}

sub _new_interval {
    my $self = shift;
    my ($bigwig,$chrom,$options) = @_;
    my $inner_iterator = $self->{inner_i};
    my %options = (%$options,
		   -seq_id => $chrom->name,
		   -start  => 1,
		   -end    => $chrom->size);
    return $inner_iterator->new($bigwig,\%options);
}

1;
