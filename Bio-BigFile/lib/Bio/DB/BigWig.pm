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

    $options{-type} ||= 'region';

    if (ref $options{-type} && ref $options{-type} eq 'ARRAY') {
	warn "This module only supports fetching one feature type at a time. Picking first one.";
	$options{-type} = $options{-type}->[0];
    }

    my $iterator_class;

    # type eq 'bin' returns a list of bins across the region. 
    # Score is ExtendedSummary hash data.
    # Start and end are beginning and end of bins
    # Select number of bins using format bin:1000,
    # defaults to 1024 bins
    if ($options{-type} =~ /^bin/) {
	$iterator_class = 'Bio::DB::BigWig::BinIterator';
    }

    # type eq 'summary' returns a feature object with a method named
    # statistical_summary(), which returns an array of extended hashes across
    # region. Start and end coordinates match request

    elsif ($options{-type} =~ /^summary/) {
	$iterator_class = 'Bio::DB::BigWig::SummaryIterator';
    }

    elsif ($options{-type} =~ /^region/) {
	$iterator_class = 'Bio::DB::BigWig::IntervalIterator';
    }

    else {
	$iterator_class = 'Bio::DB::BigWig::EmptyIterator';
    }
    
    # first deal with the problem of the user not specifying the chromosome
    return Bio::DB::BigWig::GlobalIterator->new($self,$iterator_class,\%options)
	unless $options{-seq_id};

    # now deal with the problem of the user not specifying either the
    # start or the end position
    $options{-start} ||= 1;   # that was easy!
    $options{-end}   ||= $self->bw->chromSize($options{-seq_id});
    
    return unless $options{-seq_id} && $options{-start} && $options{-end};

    return $iterator_class->new($self,\%options);
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
						 -bw     => $bw)
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

    my $bw = $self->{bw} or return;
    return $bw->bigWigSummaryArrayExtended($self->seq_id,
					   $self->start-1,
					   $self->end,
					   $bins);
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
