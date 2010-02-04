#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use Bio::Root::IO;
use FindBin '$Bin';
use constant TEST_COUNT => 19;

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
use File::Temp 'tempfile';

# low level tests (defined in lib/Bio/DB/Sam.xs) 
ok('loaded ok');

# constants
ok(defined &bbiSumMean);

my $testfile = "$Bin/../ExampleData/dpy-27-variable.bw";
my $wig      = Bio::DB::BigFile->bigWigFileOpen($testfile);
ok($wig);
ok($wig->isa('Bio::DB::BBIFile'));


# does dumping work? 
my $fh = tempfile();
ok($wig->bigWigIntervalDump('I',1,5000,0,$fh));
seek ($fh,0,0);
my $count = 0;
$count++ while <$fh>;
ok($count>1);
undef $fh;

my $interval_head = $wig->bigWigIntervalQuery('I',1,5000);
ok($interval_head);
my $nodes = 0;
for (my $i=$interval_head->head;$i;$i=$i->next) {
    $nodes++;
}
ok($nodes,$count-1);
undef $interval_head;

my $summary = $wig->bigWigSummaryArray('I',1,5000000,bbiSumMean,500);
ok($summary);
ok(ref $summary,'ARRAY');
ok(scalar @$summary,500);

my $should_fail = $wig->bigWigSummaryArray('chromFoo',1,5000000,bbiSumMean,500);
ok(!$should_fail);

my $extended_summary = $wig->bigWigSummaryArrayExtended('I',1,5000000,500);
ok ($extended_summary);
ok ($extended_summary->isa('Bio::DB::BigWigExtendedSummary'));
ok ($extended_summary->validCount(0) > 0);
ok ($extended_summary->sumData(1)/$extended_summary->validCount(1),$summary->[1]);

my $bin_sum = $wig->binStats('I',1,5000000,500);
ok ($bin_sum);
ok (scalar @$bin_sum,500);
ok($bin_sum->[300]->sumData/$bin_sum->[300]->validCount,$summary->[300]);

1;

