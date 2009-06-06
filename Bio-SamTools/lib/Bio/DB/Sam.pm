package Bio::DB::Sam;
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
    my $class        = shift;
    my %args         = @_;
    my $fa_path      = $args{-fasta} or croak "-fasta argument required";
    my $bam_path     = $args{-bam}   or croak "-bam argument required";
    my $expand_flags = $args{-expand_flags};
    -e $fa_path  && -r _  or croak "$fa_path does not exist or is not readable";
    -e $bam_path && -r _  or croak "$fa_path does not exist or is not readable";
    my $fai = Bio::DB::Sam::Fai->open($fa_path)  or croak "$fa_path open: $!";
    my $bam = Bio::DB::Bam->open($bam_path)      or croak "$bam_path open: $!";
    my $self =  bless {
	fai          => $fai,
	bam          => $bam,
	bam_path     => $bam_path,
	expand_flags => $expand_flags,
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

sub reset_read {
    my $self = shift;
    $self->{bam}->header;
}

sub targets {
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

sub ids {
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

sub maybe_build_index {
    my $self = shift;
    my $bp   = $self->{bam_path};
    unless (-e "${bp}.bai") {
	print STDERR "[bam_index_build] creating index for $bp\n";
	Bio::DB::Bam->index_build($bp);
    }
    return -e "${bp}.bai";
}

sub fetch {
    my $self     = shift;
    my $region   = shift;
    my $callback = shift;

    my $header              = $self->{bam}->header;
    my ($seqid,$start,$end) = $header->parse_region($region);
    return unless defined $seqid;
    my $index  = $self->bam_index;
    $index->fetch($self->{bam},$seqid,$start,$end,$callback,$self);
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
    my $iterator  = $args{-iterator};
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
 	return map {$self->segment($_)} $self->ids;
     }

    my %seenit;
    my @types = grep {!$seenit{$_}++} ref $types ? @$types : $types;
    @types    = 'match' unless @types;

    # the filter is intended to be inserted into a closure
    # it will return undef from the closure unless the filter
    # criteria are satisfied
    my $filter = '';
    $filter   .= $self->_filter_by_name(lc $args{-name})
	         if defined $args{-name};
    $filter   .= $self->_filter_by_attribute($args{-attributes})
 	         if defined $args{-attributes};

    # special case -- can iterate through the database using read1()
    if (@types == 1 && $types[0] =~ /^match/
	&&
	$iterator && !$use_index) {
	$self->reset_read;
	my $code = eval "sub {my \$a=shift;$filter;1}";
	die $@ if $@;
	return Bio::DB::Bam::ReadIterator->new($self->{bam},$code);
    }

    # otherwise we're going to do a little magic
    my ($features,@result);

    for my $t (@types) {

	if ($t =~ /^(match|read_pair)/) {
	    
	    # fetch the features if type is 'match' or 'read_pair'
	    $features = $self->_features($seqid,$start,$end,$filter);

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

sub _features {
    my $self = shift;
    my ($seqid,$start,$end,$filter) = @_;

    my @result;

    my $callback = defined($seqid) ? <<INDEXED : <<NONINDEXED;
sub {
    my \$a = shift;
    $filter
    return unless defined \$a->start;
    push \@result,Bio::DB::Bam::AlignWrapper->new(\$a,\$self);
    return 1;
}
INDEXED
sub {
    my \$a    = shift;
    $filter
    return 1;
}
NONINDEXED
    ;

    my $code = eval $callback;
    die $@ if $@;

    if (defined $seqid) {
 	my $region = $seqid;
 	if (defined $start) { 
 	    $region   .= ":$start";
 	    $region   .= "-$end"   if defined $end;
 	}
 	$self->fetch($region,$code);
    } 

    else {
	$self->reset_read;
	while (my $b = $self->{bam}->read1) {
 	    push @result,Bio::DB::Bam::AlignWrapper->new($b,$self) if $code->($b);
 	}
    }
    
    return \@result;
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

    return Bio::SeqFeature::Lite->new(
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
    $self->maybe_build_index;
    return Bio::DB::Bam->index_open($self->{bam_path});
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


1;
__END__
