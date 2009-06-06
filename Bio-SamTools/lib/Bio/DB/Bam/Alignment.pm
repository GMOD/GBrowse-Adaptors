package Bio::DB::Bam::Alignment;

use strict;
use warnings;
use Bio::DB::Bam::Query;
use Bio::DB::Bam::Target;
use Bio::DB::Sam::Constants;

sub get_tag_values {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    # we always expand flags
    if (my $mask = RFLAGS()->{uc $tag}) {
        # to avoid warnings when making numeric comps
	return ($self->flag & $mask) == 0 ? 0 : 1; 
    } else {
	$self->aux_get($tag);
    }
}

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

# in SAM format, alignment is always to the forward strand
sub strand {
    my $self     = shift;
    return 1;
}

sub abs_strand { shift->strand }

sub mstrand {
    my $self     = shift;
    return $self->mreversed ? -1 : 1;
}

sub display_name {
    return shift->qname;
}

sub primary_id {
    my $self = shift;
    return join ':',$self->display_name,$self->seq_id,$self->start,$self->end,$self->strand;
}

sub cigar_str {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my $result = '';
    for my $c (@$cigar) {
	my $op     = $c & BAM_CIGAR_MASK;
	my $l      = $c >> BAM_CIGAR_SHIFT();
	my $symbol = CIGAR_SYMBOLS()->[$op];
	$result .= "${symbol}${l}";
    }
    return $result;
}

sub cigar_array {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my @result;
    for my $c (@$cigar) {
	my $op     = $c & BAM_CIGAR_MASK();
	my $l      = $c >> BAM_CIGAR_SHIFT();
	my $symbol = CIGAR_SYMBOLS()->[$op];
	push @result,[$symbol,$l];
    }
    return \@result;

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

sub query {
    my $self = shift;
    return Bio::DB::Bam::Query->new($self);
}

# Target is the same as Query, but with meaning of start() and end() reversed
# for compatibility with Bio::DB::GFF and its ilk. Please use Query if you can!
sub target {
    my $self = shift;
    return Bio::DB::Bam::Target->new($self);
}

sub primary_tag { 'match'   }
sub source_tag  { 'sam/bam' }

sub hit { shift->query(@_); }


1;
