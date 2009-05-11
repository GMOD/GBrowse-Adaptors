package Bio::DB::Bam::FetchIterator;
use strict;

sub new {
    my $self = shift;
    my $list = shift;
    return bless $list,ref $self || $self;
}

sub next_seq {
    my $self = shift;
    return shift @$self;
}

1;
