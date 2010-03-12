package Bio::DB::BigWig;

#$Id$

use strict;
use warnings;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::DB::BigFile::Iterators;
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

sub bf { shift->{bw} }
sub bw { shift->{bw} }
sub fa { shift->{fa} }

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

    $start ||= 1;
    $end   ||= $self->bw->chromSize($seqid);

    return Bio::DB::BigFile::Segment->new(-bf    => $self,
					  -seq_id=> $seqid,
					  -start => $start,
					  -end   => $end);
}

sub bigwig { shift->bw }

sub get_features_by_name  { return; } # this doesn't do anything
sub get_features_by_alias { return; } # this doesn't do anything
sub get_features_by_attribute { return; } # this doesn't do anything
sub get_features_by_type {
    my $self = shift;
    return $self->features(-type=>\@_);
}
sub get_features_by_location {
    my $self = shift;
    if ($_[0] =~ /^-/) { # named argument form
	return $self->features(@_);
    } else {
	my ($seqid,$start,$end) = @_;
	return $self->features(-seq_id => $seqid,
			       -start  => $start,
			       -end    => $end);
    }

}
# kind of a cheat because there are no real IDs in bed/bigwig files
# we simply encode the location and other identifying information
# into the id
sub get_feature_by_id {
    my $self = shift;
    my $id   = shift;
    my ($chr,$start,$end) = split ':',$id;
    my @f = $self->get_features_by_location(-seq_id=>$chr,
					    -start => $start,
					    -end   => $end);
    @f or return;
    return $f[0] if @f == 1; # yay!
    
    @f = grep {$_->start == $start && $_->end == $end} @f;
    return $f[0] if @f == 1;

    warn "Did not find a single feature with id $id. Returning first one.";
    return $f[0]; # lie
}


sub seq_ids {
    my $self = shift;
    my $bw   = $self->bw;
    my $chrom_list = $bw->chromList;
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
    my %options;

    if (@_ && $_[0] !~ /^-/) {
	%options = (-type => $_[0]);
    } else {
	%options = @_;
    }

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
    my %options;

    if (@_ && $_[0] !~ /^-/) {
	%options = (-type => $_[0]);
    } else {
	%options = @_;
    }

    $options{-type} ||= 'region';

    if (ref $options{-type} && ref $options{-type} eq 'ARRAY') {
	warn "This module only supports fetching one feature type at a time. Picking first one."
	    if @{$options{-type}}>1;
	$options{-type} = $options{-type}->[0];
    }

    my $iterator_class = $self->_type_to_iterator($options{-type});
    
    # first deal with the problem of the user not specifying the chromosome
    return Bio::DB::BigFile::GlobalIterator->new($self,$iterator_class,\%options)
	unless $options{-seq_id};

    # now deal with the problem of the user not specifying either the
    # start or the end position
    $options{-start} ||= 1;   # that was easy!
    $options{-end}   ||= $self->bw->chromSize($options{-seq_id});
    
    return unless $options{-seq_id} && $options{-start} && $options{-end};

    return $iterator_class->new($self,\%options);
}

sub _type_to_iterator {
    my $self = shift;
    my $type = shift;

    return 'Bio::DB::BigWig::BinIterator'      if $type =~ /^bin/;
    return 'Bio::DB::BigWig::SummaryIterator'  if $type =~ /^summary/;
    return 'Bio::DB::BigWig::IntervalIterator' if $type =~ /^region/;
    return 'Bio::DB::BigFile::EmptyIterator';
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

use base 'Bio::DB::BigFile::IntervalIterator';

sub _query_method { 'bigWigIntervalQuery' }

sub _feature_method {'Bio::DB::BigWig::Feature'}

############################################################

package Bio::DB::BigWig::BinIterator;
use base 'Bio::DB::BigFile::BinIterator';
use Carp 'croak';

sub _query_method   { 'bigWigSummaryArrayExtended' }
sub _feature_method { 'Bio::DB::BigWig::Feature'   }


##################################################################
package Bio::DB::BigWig::SummaryIterator;

use base 'Bio::DB::BigFile::SummaryIterator';

sub _query_class { 'Bio::DB::BigWig::Summary' }


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

sub id {
    my $self = shift;
    return join (':',$self->seq_id,$self->start,$self->end);
}


##################################################################

package Bio::DB::BigWig::Summary;
use base 'Bio::DB::BigWig::Feature';

sub new {
    my $self = shift;
    my $feat = $self->SUPER::new(@_);
    my %args = @_;
    $feat->{bf} = $args{-bf} if $args{-bf};
    return $feat;
}

sub statistical_summary {
    my $self = shift;
    my $bins = shift;
    $bins ||= 1024;

    my $bf = $self->{bf}->bf or return;
    return $bf->bigWigSummaryArrayExtended($self->seq_id,
					   $self->start-1,
					   $self->end,
					   $bins);
}

##################################################################
package Bio::DB::BigFile::Segment;
use base 'Bio::DB::BigWig::Summary';

sub features {
    my $self = shift;
    return $self->{bf}->features(-seq_id => $self->seq_id,
				 -start  => $self->start,
				 -end    => $self->end,
				 -type   => $_[0] || 'region');
}

sub get_seq_stream {
    my $self = shift;
    return $self->{bf}->get_seq_stream(-seq_id => $self->seq_id,
				       -start  => $self->start,
				       -end    => $self->end,
				       -type   => $_[0] || 'region');
}

1;
