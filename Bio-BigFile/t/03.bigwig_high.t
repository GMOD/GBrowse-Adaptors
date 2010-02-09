#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use Bio::Root::IO;
use FindBin '$Bin';
use constant TEST_COUNT => 29;

use lib "$Bin/../lib","$Bin/../blib/lib","$Bin/../blib/arch";

BEGIN {
  # to handle systems with no installed Test module
  # we include the t dir (where a copy of Test.pm is located)
  # as a fallback
  eval { require Test; };
  if( $@ ) {
    use lib 't';
  }
  use Test;
  plan test => TEST_COUNT;
}

use Bio::DB::BigWig;

# high level tests
ok('loaded ok');

my $testfile = "$Bin/../ExampleData/dpy-27-variable.bw";
my $wig      = Bio::DB::BigWig->new(-bigwig=>$testfile);
ok($wig);
ok($wig->isa('Bio::DB::BigWig'));

my $iterator = $wig->get_seq_stream;
ok ($iterator);
my $nodes = 0;
while (my $f = $iterator->next_seq) {
    print $f->seq_id,":",$f->start,'..',$f->end,' ',$f->score,"\n";
    1;
}
1;

