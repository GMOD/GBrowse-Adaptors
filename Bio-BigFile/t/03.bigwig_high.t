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
my $wig      = Bio::DB::BigWig->new(-bigwig=>$testfile,
				    -fasta=>'/var/www/gbrowse2/databases/elegans_scaffolds'
    );
ok($wig);
ok($wig->isa('Bio::DB::BigWig'));

my $iterator = $wig->get_seq_stream(-seq_id=>'I',-start=>100,-end=>1000);
ok ($iterator);
my $nodes = 0;
my $inbounds = 1;
while (my $f = $iterator->next_seq) {
    $nodes++;
    $inbounds &&= $f->seq_id eq 'I' && $f->start <= 1000 && $f->end >= 100;
}
ok($nodes,11);
ok($inbounds);

my @features = $wig->features(-seq_id=>'I',-start=>100,-end=>1000);
ok (scalar @features,$nodes);

my @big_vals  = grep {$_->score >= 0.5} @features;
my @big_vals2 = $wig->features(-seq_id=>'I',-start=>100,-end=>1000,-filter=>sub {shift->score >0.5});
ok (@big_vals,@big_vals2);

# This is probably NOT how we want it to work; instead of returning one feature for each bin,
# return a single feature that knows how to produce an array of summary values
my @bins = $wig->features(-type=>'summary',-seq_id=>'I');
ok(@bins);

1;

