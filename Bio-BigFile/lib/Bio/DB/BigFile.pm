package Bio::DB::BigFile;

# $Id$

use strict;
use warnings;

=head1 NAME

Bio::DB::BigFile -- Read BigWig & BigBed files

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use Carp 'croak';
use base qw(DynaLoader Exporter);
our $VERSION = '0.01';
our @EXPORT_OK = qw(bbiSumMean bbiSumMax bbiSumMin bbiSumCoverage bbiSumStandardDeviation);
our @EXPORT = @EXPORT_OK;

bootstrap Bio::DB::BigFile;

use constant bbiSumMean => 0;
use constant bbiSumMax  => 1;
use constant bbiSumMin  => 2;
use constant bbiSumCoverage => 3;
use constant bbiSumStandardDeviation => 4;

package Bio::DB::BBIFile;

sub binStats {
    my $self = shift;
    my $extended_summary = $self->bigWigSummaryArrayExtended(@_);
    my @tie;
    tie @tie,'Bio::DB::BigWig::binStats',$extended_summary;
    return \@tie;
}

package Bio::DB::BigWig::binStats;

use base 'Tie::Array';
use Carp 'croak';

sub TIEARRAY {
    my $class = shift;
    my $summary = shift;
    $summary 
	or croak "Usage: tie(\@array,'$class',\$summary), where \$summary is a Bio::DB::BigWigExtendedSummary object";
    return bless \$summary,ref($class) || $class;
}

sub FETCH {
    my $self  = shift;
    my $index = shift;
    return Bio::DB::BigWig::binStatElement->new($$self,$index);
#    return $$self->validCount($index);
}

sub FETCHSIZE {
    my $self = shift;
    return $$self->size;
}

package Bio::DB::BigWig::binStatElement;

sub new {
    my $class = shift;
    my ($base,$index) = @_;
    return bless [$base,$index],ref $class || $class;
}

sub validCount {
    my $self = shift;
    $self->[0]->validCount($self->[1]);
}

sub minVal {
    my $self = shift;
    $self->[0]->minVal($self->[1]);
}

sub maxVal {
    my $self = shift;
    $self->[0]->maxVal($self->[1]);
}

sub sumData {
    my $self = shift;
    $self->[0]->sumData($self->[1]);
}

sub sumSquares {
    my $self = shift;
    $self->[0]->sumSquares($self->[1]);
}

=head1 METHODS

=head1 EXAMPLES

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::Graphics>, L<Bio::Graphics::Browser2>

=head1 AUTHOR

Lincoln Stein E<lt>lincoln.stein@oicr.on.caE<gt>.
E<lt>lincoln.stein@bmail.comE<gt>

Copyright (c) 2009 Ontario Institute for Cancer Research.

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut

1;



