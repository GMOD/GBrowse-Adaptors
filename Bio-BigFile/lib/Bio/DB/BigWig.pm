package Bio::DB::BigWig;

#$Id$

use strict;
use warnings;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::Graphics::Feature;
use Carp 'croak';

# high level interface to BigWig files

sub new {
    my $self = shift;
    my %args = @_;

    my $bw_path       = $args{-bigwig}  or croak "-bigwig argument required";
    my $fa_path       = $args{-fasta};
    my $dna_accessor  = $self->new_dna_accessor($fa_path);
    
    unless ($self->is_remote($bw_path)) {
	-e $bw_path or croak "$bw_path does not exist";
	-r _  or croak "is not readable";
    }

    my $bw = Bio::DB::BigFile->bigWigFileOpen($bw_path)
	or croak "$bw_path open: $!";

    return bless {
	bw => $bw,
	fa => $dna_accessor
    },ref $self || $self;
}

sub bw { shift->{bw} }
sub fa { shift->{fa} }

sub bigwig { shift->bw }

sub seq_ids {
    my $self = shift;
    my $bw   = $self->bw;
    my $chrom_list = $bw->bbiChromInfoHead;
    my @list;
    for (my $c=$chrom_list->head;$c;$c=$c->next) {
	push @list,$c->name;
    }
    return @list;
}

sub length {
    my $self = shift;
    my $seqid = shift;
    return $self->bw->chromSize($seqid);
}

sub features {
    my $self    = shift;
    my %options = @_;

    my $iterator = $self->get_seq_stream(%options);
    return $iterator if $options{-iterator};
    
    my @result;
    while (my $f = $iterator->next_seq) {
	push @result,$f;
    }

    return @result;
}

sub get_seq_stream {
    my $self    = shift;
    my %options = @_;

    return $self->_get_summary_stream(\%options)
	if $options{-type} && $options{-type} =~ /summary/;
    
    # first deal with the problem of the user not specifying the chromosome
    return Bio::DB::BigWig::GlobalIterator->new($self,\%options)
	unless $options{-seq_id};

    # now deal with the problem of the user not specifying either the
    # start or the end position
    $options{-start} ||= 1;   # that was easy!
    $options{-end}   ||= $self->bw->chromSize($options{-seq_id});

    return unless $options{-seq_id} && $options{-start} && $options{-end};

    return Bio::DB::BigWig::IntervalIterator->new($self,\%options);
    
}

sub _get_summary_stream {
    my $self = shift;
    my $options = shift;
    my ($type,$bins) = $options->{-type} =~ /^(summary)(?::(\d+))?/i
	or croak "invalid call to _get_summary_stream. -type argument must be summary[:bins]";

    $bins ||= 1024;

    # first deal with the problem of the user not specifying the chromosome
    return Bio::DB::BigWig::GlobalSummaryIterator->new($self,$bins,$options)
	unless $options->{-seq_id};

    # now deal with the problem of the user not specifying either the
    # start or the end position
    $options->{-start} ||= 1;   # that was easy!
    $options->{-end}   ||= $self->bw->chromSize($options->{-seq_id});

    return unless $options->{-seq_id} && $options->{-start} && $options->{-end};

    return Bio::DB::BigWig::SummaryIterator->new($self,$bins,$options);
}

sub seq {
    my $self = shift;
    my ($seqid,$start,$end) = @_;
    my $fa   = $self->fa;
    return $fa ? $fa->seq($seqid,$start,$end) : 'N' x ($end-$start+1);
}

sub is_remote {
    my $self = shift;
    my $path = shift;
    return $path =~ /^(http|ftp):/;
}

sub new_dna_accessor {
    my $self     = shift;
    my $accessor = shift;

    return unless $accessor;

    if (-e $accessor) {  # a file, assume it is a fasta file
	eval "require Bio::DB::Fasta" unless Bio::DB::Fasta->can('new');
	my $a = Bio::DB::Fasta->new($accessor)
	    or croak "Can't open FASTA file $accessor: $!";
	return $a;
    }

    if (ref $accessor && UNIVERSAL::can($accessor,'seq')) {
	return $accessor;  # already built
    }

    my $obj = eval $accessor;  # maybe we got some code????
    if ($obj && ref $obj && UNIVERSAL::can($obj,'seq')) {
	return $obj;
    }

    return;
}

############################################################

package Bio::DB::BigWig::IntervalIterator;

sub new {
    my $self   = shift;
    my ($bigwig,$options) = @_;
    my $bw     = $bigwig->bw;
    my $head   = $bw->bigWigIntervalQuery($options->{-seq_id},
					  $options->{-start}-1,
					  $options->{-end})
	or return;
    return bless {
	head    => $head,   # keep in scope so not garbage collected
	seq_id  => $options->{-seq_id},
	current => $head->head,
	bigwig  => $bigwig,
	options => $options,
    },ref $self || $self;
}

sub next_seq {
    my $self = shift;
    my $filter = $self->{options}{-filter};

    my ($i,$f);

    for ($i = $self->{current};$i;$i=$i->next) {
	$f = Bio::DB::BigWig::Feature->new(-seq_id => $self->{seq_id},
					   -start  => $i->start+1,
					   -end    => $i->end,
					   -score  => $i->value,
					   -type   => 'region',
					   -fa     => $self->{bigwig}->fa,
	    );
	last if !$filter || $filter->($f);
    }

    if ($i) {
	$self->{current} = $i->next;
	return $f;
    }
    else {
	$self->{current} = undef;
	return;
    }
}

############################################################

package Bio::DB::BigWig::SummaryIterator;

sub new {
    my $self   = shift;
    my ($bigwig,$bins,$options) = @_;
    my $bw     = $bigwig->bw;
    my $arry   = $bw->bigWigSummaryArrayExtended($options->{-seq_id},
						 $options->{-start}-1,
						 $options->{-end},
						 $bins)
	or return;
    my $binsize = ($options->{-end}-$options->{-start}+1)/$bins;
    return bless {
	array   => $arry,
	bigwig  => $bigwig,
	start   => $options->{-start},
	binsize => $binsize,
    },ref $self || $self;
}

sub next_seq {
    my $self = shift;
    my $filter = $self->{options}{-filter};

    my $array  = $self->{array};
    my $i      = shift @$array;
    my $f;
    while ($i) {
	$f = Bio::DB::BigWig::Feature->new(-seq_id => $self->{seq_id},
					   -start  => int($self->{start}),
					   -end    => int($self->{start}+$self->{binsize}),
					   -score  => $i,
					   -type   => 'summary',
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
    my $fa     = $self->{fa} or return;
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

package Bio::DB::BigWig::GlobalIterator;

sub new {
    my $self = shift;
    my ($bigwig,$options) = @_;
    my $cl = $bigwig->bw->chromList or return;
    return bless {
	cl_head  => $cl,     # keep in scope so not garbage collected
	current  => $cl->head,
	bigwig   => $bigwig,
	options  => $options,
	interval => $self->_new_interval($bigwig,$cl->head,$options),
    },ref $self || $self;
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
    my %options = (%$options,
		   -seq_id => $chrom->name,
		   -start  => 0,
		   -end    => $chrom->size);
    return Bio::DB::BigWig::IntervalIterator->new($bigwig,\%options);
}


##################################################################

package Bio::DB::BigWig::GlobalSummaryIterator;

sub new {
    my $self = shift;
    my ($bigwig,$bins,$options) = @_;
    my $cl = $bigwig->bw->chromList or return;
    return bless {
	cl_head  => $cl,     # keep in scope so not garbage collected
	current  => $cl->head,
	bigwig   => $bigwig,
	options  => $options,
	bins     => $bins,
	interval => $self->_new_summary_array($bigwig,$cl->head,$bins,$options),
    },ref $self || $self;
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
	$self->{interval} = $self->_new_summary_array($self->{bigwig},
						      $self->{current},
						      $self->{bins},
						      $self->{options});
	my $next = $self->{interval}->next_seq;
	return $next if $next;
    }
}

sub _new_summary_array {
    my $self = shift;
    my ($bigwig,$chrom,$options) = @_;
    my %options = (%$options,
		   -seq_id => $chrom->name,
		   -start  => 0,
		   -end    => $chrom->size);
    return Bio::DB::BigWig::SummaryIterator->new($bigwig,\%options);
}


1;
