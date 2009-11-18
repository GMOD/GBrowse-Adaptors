package Bio::DB::Bam::Target;
use strict;

use base 'Bio::DB::Bam::Query';

sub start {
    my $self = shift;
    return $self->strand > 0 ? $self->low : $self->high;
}

sub end {
    my $self = shift;
    return $self->strand > 0 ? $self->high : $self->low;
}


1;

