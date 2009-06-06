package Bio::DB::Bam::AlignWrapper;

use strict;

*RFLAGS = \&Bio::DB::Bam::Alignment::RFLAGS;

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
  my $self = shift or die "autoload called for non-object symbol $func_name";
  $self->{align}->$func_name(@_);
}

sub can {
    my $self = shift;
    return 1 if $self->SUPER::can(@_);
    return $self->{align}->can(@_);
}
sub seq_id {
    my $self = shift;
    my $tid  = $self->tid;
    $self->{sam}->target_name($tid);
}

sub expand_flags {
    shift->{sam}->expand_flags(@_);
}

sub primary_tag { return 'match' }
sub primary_id {
    my $self = shift;
    return join ';',
       map {s/;/%3B/g; $_}
       ($self->display_name,
	$self->seq_id,
	$self->start,
	$self->end,
	$self->strand);
}
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

sub seq      {
    my $self   = shift;
    return Bio::PrimarySeq->new(-seq => $self->dna,
				-id  => $self->seq_id);
}

sub subseq {
    my $self = shift;
    my ($start,$end) = @_;
    $start = 1 if $start < 1;
    $end   = $self->high if $end > $self->high;
    return Bio::PrimarySeq->new(-seq=>substr($self->dna,
					     $start-1,
					     $end-$start+1)
				);
}

sub dna {
    my $self = shift;
    my $region = $self->seq_id.':'.$self->start.'-'.$self->end;
    return $self->{sam}->fai->fetch($region);
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
    my @flag_tags = $self->expand_flags ? keys %{RFLAGS()} : 'FLAGS';
    return (@aux_tags,@flag_tags);
}

sub get_tag_values {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if ($self->expand_flags && 
	(my $mask = RFLAGS()->{uc $tag})) {  # special tag
        # to avoid warnings when making numeric comps
	return ($self->flag & $mask) == 0 ? 0 : 1; 
    } elsif ($tag eq 'FLAGS') {
	$self->flag_str;
    } else {
	$self->aux_get($tag);
    }
}

sub has_tag {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if ($self->expand_flags && 
	(my $mask = RFLAGS()->{uc $tag})) {  # special tag
	return 1;
    } elsif ($tag eq 'FLAGS') {
	return 1;
    } else {
	my %keys = map {$_=>1} $self->aux_keys;
	return exists $keys{uc $tag};
    }
}

1;
