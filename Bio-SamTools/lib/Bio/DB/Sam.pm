package Bio::DB::Sam;
use strict;
use warnings;
use Carp 'croak';
use Bio::SeqFeature::Lite;

use base 'DynaLoader';
our $VERSION = '0.01';

bootstrap Bio::DB::Sam;

sub new {
    my $class     = shift;
    my %args      = @_;
    my $fa_path   = $args{-fasta} or croak "-fasta argument required";
    my $bam_path  = $args{-bam}   or croak "-bam argument required";
    -e $fa_path  && -r _  or croak "$fa_path does not exist or is not readable";
    -e $bam_path && -r _  or croak "$fa_path does not exist or is not readable";
    my $fai = Bio::DB::Sam::Fai->open($fa_path)  or croak "$fa_path open: $!";
    my $bam = Bio::DB::Bam->open($bam_path)      or croak "$bam_path open: $!";
    my $self =  bless {
	fai      => $fai,
	bam      => $bam,
	bam_path => $bam_path
    },ref $class || $class;
    $self->header;  # catch it
    return $self;
}

sub header {
    my $self = shift;
    return $self->{header} ||= $self->{bam}->header;
}

sub fai { shift->{fai} }

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
    $types        = [$args{-class}] if !@$types && defined $args{-class};
    my $use_index = defined $seqid;

    # we do some special casing to retrieve target (reference) sequences
    # if they are requested
    if (defined $args{-name} 
	&& (!@$types || "@$types" =~ /region|chromosome/)) {
	my @results = $self->_segment_search(lc $args{-name});
	warn "results = @results";
	return @results if @results;
    } elsif ($types->[0] =~ /region|chromosome/) {
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
	return Bio::DB::Sam::ReadIterator->new($self->{bam},$code);
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

    return $iterator ? Bio::DB::Sam::FetchIterator->new(\@result)
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
    return 0 unless defined \$a->start;
    push \@result,Bio::DB::Bam::AlignWrapper->new(\$a,\$self);
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

    if (my $regexp = $self->_glob_match($name)) {
	return "return unless \$a->qname =~ /^$regexp\$/i;\n";
    } else {
	return "return unless lc \$a->qname eq '$name';\n";
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

package Bio::DB::Bam::AlignWrapper;

our $AUTOLOAD;

sub new {
    my $package = shift;
    my ($align,$sam) = @_;
    return bless {sam   => $sam,
		  align => $align},ref $package || $package;
}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  return if $func_name eq 'DESTROY';
  my $self = shift or die;
  $self->{align}->$func_name(@_);
}

sub can {
    shift->{align}->can(@_);
}
sub seq_id {
    my $self = shift;
    my $tid  = $self->tid;
    $self->{sam}->target_name($tid);
}

sub primary_tag { return 'match' }
sub abs_ref    { shift->seq_id }
sub abs_start  { shift->start  }
sub abs_end    { shift->end    }
sub low        { shift->start  }
sub high       { shift->end    }
sub type       { shift->primary_tag }
sub method     { shift->primary_tag }
sub source_tag { return 'sam/bam'; }
sub source     { return shift->source_tag; }
sub name       { shift->qname }
sub class      { shift->primary_tag }

# required by API
sub get_SeqFeatures { return }

package Bio::DB::Bam::Alignment;

use constant CIGAR_SYMBOLS   => [qw(M I D N S H P)];
use constant BAM_CIGAR_SHIFT => 4;
use constant BAM_CIGAR_MASK  => (1 << BAM_CIGAR_SHIFT) - 1;
use constant BAM_CMATCH      => 0;
use constant BAM_CINS        => 1;
use constant BAM_CDEL        => 2;
use constant BAM_CREF_SKIP   => 3;
use constant BAM_CSOFT_CLIP  => 4;
use constant BAM_CHARD_CLIP  => 5;
use constant BAM_CPAD        => 6;

use constant FLAGS => {
    0x0001 => 'PAIRED',
    0x0002 => 'MAP_PAIR',
    0x0004 => 'UNMAPPED',
    0x0008 => 'M_UNMAPPED',
    0x0010 => 'REVERSED',
    0x0020 => 'M_REVERSED',
    0x0040 => 'FIRST_MATE',
    0x0080 => 'SECOND_MATE',
    0x0100 => 'NOT_PRIMARY',
    0x0200 => 'QC_FAILED',
    0x0400 => 'DUPLICATE'
};

use constant RFLAGS => {reverse %{FLAGS()}};

sub start {
    my $self = shift;
    return if $self->unmapped;
    return $self->pos+1;
}

sub end {
    my $self = shift;
    return if $self->unmapped;
    return $self->calend;
}

sub stop { shift->end }

sub strand {
    my $self     = shift;
    return $self->reversed ? -1 : 1;
}

sub mstrand {
    my $self     = shift;
    return $self->mreversed ? -1 : 1;
}

sub display_name {
    return shift->qname;
}

sub attributes {
    my $self = shift;
    my $tag  = shift;
    if (defined $tag) {
	return $self->get_tag_values($tag);
    } else {
	return map {$_=>$self->get_tag_values($_)} $self->get_all_tags;
    }
}

sub get_all_tags {
    my $self      = shift;
    my @aux_tags  = $self->aux_keys;
    my @flag_tags = keys %{RFLAGS()};
    return (@aux_tags,@flag_tags);
}

sub get_tag_values {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if (my $mask = RFLAGS()->{uc $tag}) {  # special tag
        # to avoid warnings when making numeric comps
	return ($self->flag & $mask) == 0 ? 0 : 1; 
    } else {
	$self->aux_get($tag);
    }
}

sub has_tag {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if (my $mask = RFLAGS()->{uc $tag}) {  # special tag
	return 1;
    } else {
	my %keys = map {$_=>1} $self->aux_keys;
	return exists $keys{uc $tag};
    }
}

sub cigar_str {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my $result = '';
    for my $c (@$cigar) {
	my $op     = $c & BAM_CIGAR_MASK;
	my $l      = $c >> BAM_CIGAR_SHIFT;
	my $symbol = CIGAR_SYMBOLS->[$op];
	$result .= "${symbol}${l}";
    }
    return $result;
}

sub flag_str {
    my $self  = shift;
    my $flag  = $self->flag;
    my $flags = FLAGS;
    return join '|',map {$flags->{$_}}
			 grep {$flag & $_}
			 sort {$a<=>$b}
			 keys %{$flags};
}

sub length {
    my $self = shift;
    return $self->end-$self->start+1;
}

sub mate_start {
    shift->mpos+1;
}

sub mate_len {
    my $self    = shift;
    my $ins_len = $self->isize or return;
    my $len     = $self->length;

    my $adjust = 0;
    my @cigar   = $self->cigar_str =~ /(\w)(\d+)/g;
    while (@cigar) {
	my ($op,$len) = splice(@cigar,0,2);
	$adjust += $len if $op eq 'I';
	$adjust -= $len if $op eq 'D';
    }

    return $adjust + $ins_len + ($self->start-$self->mate_start) if $ins_len > 0;
    return $adjust + $self->mate_start-($self->start+$ins_len)   if $ins_len < 0;
}

sub mate_end {
    my $self = shift;
    return unless $self->mate_len;
    return $self->mate_start+$self->mate_len-1;
}

package Bio::DB::Sam::ReadIterator;
sub new {
    my $self = shift;
    my ($bam,$filter) = @_;
    return bless {bam=>$bam,
		  filter=>$filter},ref $self || $self;
}
sub next_seq {
    my $self = shift;
    while (my $b = $self->{bam}->read1) {
	return $b if $self->{filter}->($b);
    }
    return;
}

package Bio::DB::Sam::FetchIterator;

sub new {
    my $self = shift;
    my $list = shift;
    return bless $list,ref $self || $self;
}

sub next_seq {
    my $self = shift;
    return shift @$self;
}

package Bio::DB::Sam::Segment;

use Bio::PrimarySeq;

sub new {
    my $class                   = shift;
    my ($db,$seqid,$start,$end) = @_;
    return bless {
	db     => $db,
	seqid  => $seqid,
	start  => $start,
	end    => $end},ref $class || $class;
}

sub db       { shift->{db}    };

# required by api
sub seq_id   { shift->{seqid} };
# required by api
sub start    { shift->{start} };
# required by api
sub end      { shift->{end}   };
# required by api
sub strand   { 0              };
# required by api
sub length   { 
    my $self = shift;
    return $self->end - $self->start + 1;
}
# required by api
sub seq      {
    my $self   = shift;
    my $db     = $self->db;
    my $region = $self->seq_id.':'.$self->start.'-'.$self->end;
    return Bio::PrimarySeq->new(-seq => $db->fai->fetch($region),
				-id  => $self->seq_id);
}
# required by api
sub primary_tag {
    my $self = shift;
    return 'region';
}
# required by api
sub source_tag { return 'sam/bam' }
# required by api
sub name    { shift->seq_id }
# required by api
sub factory { shift->db  }
# required by api
sub display_name { shift->name  }
# required by api
sub get_SeqFeatures { return;   }
# required by api
sub method { shift->primary_tag }
# required by api
sub get_tag_values {  return; }
# required by api
sub score { return;  }
# required by api
sub class { 'sequence'  }
# required by api
sub abs_ref   { shift->seq_id }
sub abs_start { shift->start  }
sub abs_end   { shift->end    }


1;
__END__
