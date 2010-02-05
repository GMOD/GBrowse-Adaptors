#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use Bio::Root::IO;
use FindBin '$Bin';
use constant TEST_COUNT => 12;

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

use Bio::DB::BigFile;

my $bed = Bio::DB::BigFile->bigBedFileOpen("$Bin/../ExampleData/refSeqTest.flat.bb");
ok($bed);

my $interval_list = $bed->bigBedInterval('chr1',0=>20_000_000);
ok($interval_list);
ok(ref $interval_list->head,'Bio::DB::BigBedInterval');
my $nodes = 0;
for (my $i=$interval_list->head;$i;$i=$i->next) {
    $nodes++;
}
ok($nodes,6);

my $summary = $bed->bigBedSummaryArray('chr1',0,12_000_000,bbiSumCoverage,500);
ok($summary);
ok(ref $summary,'ARRAY');
ok(scalar @$summary,500);
ok(!defined $summary->[0]);
ok(defined $summary->[-1]);

my $es = $bed->bigBedSummaryArrayExtendedHash('chr1',0,12_000_000,500);
ok ($es);
ok(scalar @$es,500);
ok($es->[-2]{validCount}/(12_000_000/500),$summary->[-2]);


1;

