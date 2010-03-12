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
    
    my $bb_path       = $args{-bigbed}  or croak "-bigbed argument required";
    my $fa_path       = $args{-fasta};
    my $dna_accessor  = $self->new_dna_accessor($fa_path);
    
    unless ($self->is_remote($bb_path)) {
	-e $bb_path or croak "$bb_path does not exist";
	-r _  or croak "is not readable";
    }

    my $bb = Bio::DB::BigFile->bigBedFileOpen($bb_path)
	or croak "$bb_path open: $!";

    return bless {
	bw => $bb,
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
    return 'Bio::DB::BigBed::FeatureIterator'  if $type =~ /^(region|feature|mRNA)/;
    return 'Bio::DB::BigFile::EmptyIterator';
}

# kind of a cheat because there are no real IDs in bed/bigwig files
# we simply encode the location and other identifying information
# into the id
sub get_feature_by_id {
    my $self = shift;
    my $id   = shift;
    my ($name,$chr,$start,$end,$strand,$parts) = split ':',$id;
    my @f = $self->get_features_by_location(-seq_id=>$chr,
					    -start => $start,
					    -end   => $end,
					    -strand=> $strand);
    @f or return;
    return $f[0] if @f == 1; # yay!

    @f = grep { $_->display_name eq $name } @f    if defined $name;
    return $f[0] if @f == 1;

    @f = grep { $_->get_SeqFeatures == $parts} @f if $parts;
    return $f[0] if @f == 1;

    warn "Did not find a single feature with id $id. Returning first one.";
    return $f[0]; # lie
}



############################################################

package Bio::DB::BigBed::FeatureIterator;
use Carp 'croak';

use base 'Bio::DB::BigFile::IntervalIterator';

sub _query_method   { 'bigBedIntervalQuery'      }

sub _feature_method { 'Bio::DB::BigBed::Feature' }

sub _make_feature {
    my $self     = shift;
    my $raw_item = shift;
    my $method   = $self->_feature_method;
    my ($name,$score,$strand,
	$thickStart,$thickEnd,$itemRGB,
	$blockCount,$blockSizes,$blockStarts)  = split /\s+/,$raw_item->rest;

    $strand = defined $strand ? $strand eq '-' ? -1 
	                                       : +1 
			      : 0;

    my $type = 'region';
    my @children;
    if ($blockCount && $blockCount > 0) {
	$type    = 'mRNA';
	my $bits = $self->split_gene_bits($raw_item->start,$raw_item->end,$strand,
					  $thickStart,$thickEnd,
					  $blockCount,$blockSizes,$blockStarts);
	my $fa = $self->{bigfile}->fa;
	@children= map {
	    $method->new(-fa     => $fa,
			 -start  => $_->[0],
			 -end    => $_->[1],
			 -type   => $_->[2],
			 -strand => $strand);
	} @$bits;
    }
    

    my @args = (-seq_id => $self->{seq_id},
		-start  => $raw_item->start+1,
		-end    => $raw_item->end,
		-type   => $type,
		-strand => $strand,
		-fa     => $self->{bigfile}->fa);
    push @args,(-display_name => $name)             if defined $name;
    push @args,(-score        => $score)            if defined $score;
    push @args,(-attributes   => {RGB => $itemRGB}) if defined $itemRGB;
    push @args,(-segments     => \@children)        if @children;
    
    return $method->new(@args);
}


sub split_gene_bits {
    my $self = shift;
    my ($chromStart,$chromEnd,$strand,
	$thickStart,$thickEnd,
	$numBlocks,$blockSizes,
	$blockStarts) = @_;

    my ($leftUTR,$rightUTR) = $strand >= 0 ? ('five_prime_UTR', 'three_prime_UTR') 
					   : ('three_prime_UTR','five_prime_UTR');

    # no internal structure, so just create UTRs and one CDS in the middle
    # remember that BED format uses 0-based indexing, hence the +1s
    unless ($blockSizes) {  
	my @bits = ([$chromStart+1,$thickStart,$leftUTR],
		    [$thickStart+1,$thickEnd,'CDS'],
		    [$thickEnd+1,$chromEnd,$rightUTR]);
	return \@bits;
    }

    # harder -- we have internal exons
    my @block_sizes  = split ',',$blockSizes;
    my @block_starts = split ',',$blockStarts;
    croak "Invalid BED file: blockSizes != blockStarts"
	unless @block_sizes == @block_starts && @block_sizes == $numBlocks;

    my @bits;
    for (my $i=0;$i<@block_starts;$i++) {
	my $start = $chromStart + $block_starts[$i];	
	my $end   = $chromStart + $block_starts[$i] + $block_sizes[$i];

	if ($start < $thickStart) {
	    if ($end < $thickStart) {          # UTR wholly contained in an exon
		push @bits,[$start+1,$end,$leftUTR];
	    }
	    elsif ($end >= $thickStart) {      # UTR partially contained in an exon
		push @bits,[$start+1,$thickStart,$leftUTR];
		push @bits,[$thickStart+1,$end,'CDS'];
	    }
	}

	elsif ($start < $thickEnd) {
	    if ($end <= $thickEnd) {           # CDS wholly contained in an exon
		push @bits,[$start+1,$end,'CDS'];
	    }
	    elsif ($end > $thickEnd) {         # CDS partially contained in an exon
		push @bits,[$start+1,$thickEnd,'CDS'];
		push @bits,[$thickEnd+1,$end,$rightUTR];
	    }
	}

	elsif ($start > $thickEnd) {
	    push @bits,[$start+1,$end,$rightUTR];  # UTR wholly contained in an exon
	}

	else {
	    croak "Programmer error when calculating UTR bounds";
	}

    }

    return \@bits;
}

############################################################

package Bio::DB::BigBed::Feature;
use base 'Bio::DB::BigWig::Feature';

sub id {
    my $self = shift;
    my $chr    = $self->seq_id;
    my $start  = $self->start;
    my $end    = $self->end;
    my $strand = $self->strand;
    my $parts  = $self->get_SeqFeatures;
    my $name   = $self->display_name;
    return join ':',$name,$chr,$start,$end,$strand,$parts;
}

############################################################

package Bio::DB::BigBed::BinIterator;
use base 'Bio::DB::BigFile::BinIterator';

sub _query_method   { 'bigBedSummaryArrayExtended' }
sub _feature_method {'Bio::DB::BigBed::Feature'    }


##################################################################
package Bio::DB::BigBed::SummaryIterator;
use base 'Bio::DB::BigFile::SummaryIterator';

sub _query_class { 'Bio::DB::BigBed::Summary' }


##################################################################
package Bio::DB::BigBed::Summary;
use base 'Bio::DB::BigWig::Summary';

sub statistical_summary {
    my $self = shift;
    my $bins = shift;
    $bins ||= 1024;

    my $bf = $self->{bf}->bf or return;
    return $bf->bigBedSummaryArrayExtended($self->seq_id,
					   $self->start-1,
					   $self->end,
					   $bins);
}

1;
