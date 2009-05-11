package Bio::DB::Bam::ReadIterator;

use strict;

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

1;
