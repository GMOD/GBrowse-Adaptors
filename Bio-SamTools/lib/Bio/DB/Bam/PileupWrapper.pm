package Bio::DB::Bam::PileupWrapper;
#$Id: PileupWrapper.pm,v 1.1 2009-06-23 08:30:38 lstein Exp $
use strict;
use Bio::DB::Bam::AlignWrapper;

our $AUTOLOAD;
use Carp 'croak';

sub new {
    my $package = shift;
    my ($align,$sam) = @_;
    return bless {sam    => $sam,
		  pileup => $align},ref $package || $package;

}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  return if $func_name eq 'DESTROY';

  no strict 'refs';
  $_[0] or die "autoload called for non-object symbol $func_name";
  croak qq(Can't locate object method "$func_name" via package "$pack")
      unless $_[0]->{pileup}->can($func_name);

  *{"${pack}::${func_name}"} = sub { shift->{pileup}->$func_name(@_) };

  shift->$func_name(@_);
}

sub can {
    my $self = shift;
    return 1 if $self->SUPER::can(@_);
    return $self->{pileup}->can(@_);
}

sub alignment {
    my $self = shift;
    return Bio::DB::Bam::AlignWrapper->new($self->{pileup}->b,$self->{sam});
}

1;

