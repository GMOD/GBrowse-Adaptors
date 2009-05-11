package Bio::DB::Bam::Query;

use strict;
use Bio::DB::Sam;

use constant CIGAR_SKIP      => {BAM_CREF_SKIP  => 1,
 				 BAM_CSOFT_CLIP => 1,
 				 BAM_CHARD_CLIP => 1};


sub new {
    my $self      = shift;
    my $alignment = shift;
    bless \$alignment,ref $self || $self;
}

sub seq_id {
    my $self = shift;
    $$self->qname;
}

sub subseq {
    my $self = shift;
    my ($start,$end) = @_;
    $start = 1 if $start < 1;
    $end   = $self->high if $end > $self->high;
    ($end,$start) = ($start,$end) if $start > $end;
    return Bio::PrimarySeq->new(-seq=>substr($self->dna,
					     $start-1,
					     $end-$start+1)
				);
}

sub primary_tag { ${shift()}->primary_tag }
sub source_tag  { ${shift()}->source_tag  }

sub start {
    my $self = shift;
    return $self->low;
}

sub end {
    my $self = shift;
    return $self->high;
}

sub low {
    my $self       = shift;
    my $cigar_arry = $$self->cigar_array;
    my $start      = 1;
    for my $c (@$cigar_arry) {
	last unless CIGAR_SKIP->{$c->[0]};
	$start += $c->[1];
    }
    $start;
}

sub high {
    my $self      = shift;
    my $len       = $$self->cigar2qlen;
    my $cigar_arry = $$self->cigar_array;

    # alignment stops at first non-clip CIGAR position
    my $i = $len - 1;
    for my $c (reverse @$cigar_arry) {
	last unless CIGAR_SKIP->{$c->[0]};
	$len -= $c->[1];
    }
    return $len;
}

sub length {
    my $self = shift;
    $$self->cigar2qlen;
}

sub name {
    my $self = shift;
    $$self->qname;
}

sub display_name {shift->name}

sub seq {
    my $self = shift;
    my $dna  = $self->strand > 0 ? $$self->qseq : reversec($$self->qseq);
    return Bio::PrimarySeq->new(-seq => $dna,
				-id  => $$self->qname);
}

sub dna {
    my $self = shift;
    return $$self->qseq;
}

sub strand { 
    my $self = shift;
    return $$self->reversed ? -1 : 1;
}

sub reversec {
    my $dna = shift;
    $dna =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}


1;
