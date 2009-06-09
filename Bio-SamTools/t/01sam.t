#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use Bio::Root::IO;
use FindBin '$Bin';
use constant TEST_COUNT => 57;

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

use Bio::DB::Sam;

# low level tests (defined in lib/Bio/DB/Sam.xs) 
{
    my $bamfile = "$Bin/data/ex1.bam";
    my $bam     = Bio::DB::Bam->open($bamfile);
    ok($bam);

    my $header  = $bam->header;
    my $targets = $header->n_targets;
    ok($targets,2);

    my $target_names = $header->target_name;
    ok($target_names);
    ok(scalar @$target_names,2);
    ok($target_names->[0],'seq1');
    
    my $target_lens = $header->target_len;
    ok($target_lens);
    ok(scalar @$target_lens,2);
    ok($target_lens->[0],1575);
    
    my $text = $header->text;
    ok(length $text > 0);
    
    my $fai  = Bio::DB::Sam::Fai->open("$Bin/data/ex1.fa");
    my $seq  = $fai->fetch('seq2:51-1000');
    ok(length $seq,950);

    my $count;
    while (my $b = $bam->read1) {
	$count++;
    }
    ok($count,3307);
    
    my @result = $header->parse_region('seq2:51-1000');
    ok($result[0],1);
    ok($result[1],50);
    @result    = $header->parse_region('seq_invalid:51-1000');
    ok(scalar @result,0);
    
    my $index = Bio::DB::Bam->index_open($bamfile);
    ok($index);

    my @a;
    my $print_region = sub {
	my($alignment,$data) = @_;
	push @a,$alignment;
	return;
    };

    $index->fetch($bam,$header->parse_region('seq2'),$print_region,"foobar");
    ok(scalar @a > 1);


    # try to get coverage
    my $coverage = $index->coverage($bam,
				    $header->parse_region('seq2'),
				    100);
    ok(scalar @$coverage,100);
    my @c = sort {$a<=>$b} @$coverage;
    ok($c[0]  >= 0);
    ok($c[-1] < 1000);
}

# high level tests (defined in lib/Bio/DB/Sam.pm)
{
    my $sam = Bio::DB::Sam->new(-fasta=>"$Bin/data/ex1.fa",
			        -bam  =>"$Bin/data/ex1.bam",
				-expand_flags => 1,
	);
    ok($sam);
    ok($sam->targets,2);

    ok($sam->length('seq1'),1575);
    ok(join $sam->ids,'seq1 seq2');

    my $seg = $sam->segment('seq1');
    ok($seg);
    ok($seg->length,1575);
    my $seq = $seg->seq;
    ok($seq->isa('Bio::PrimarySeq'));
    ok(length $seq->seq,1575);

    my $dummy = eval {Bio::DB::Sam->new(-fasta=>"invalid_path.txt",
					-bam  =>"invalid_path.txt")};
    ok($dummy,undef);
    ok($@ =~ /does not exist/);
    
    my @alignments = 
	$sam->get_features_by_location(
	    -seq_id => 'seq2',
	    -start  => 500,
	    -end    => 800
	);
    ok(scalar @alignments,442);
    ok($alignments[0]->seq_id,'seq2');

    ok(length scalar $alignments[0]->qscore,length $alignments[0]->dna);

    my @keys = $alignments[0]->get_all_tags;
    ok(scalar @keys,17);
    ok($alignments[0]->get_tag_values('MF'),18);

    my %att  = $alignments[0]->attributes;
    ok(scalar(keys %att),17);
    ok($alignments[0]->cigar_str,'M35');

    $sam->expand_flags(0);
    @keys = $alignments[0]->get_all_tags;
    ok(scalar @keys,7);

    my $query = $alignments[0]->query;
    ok($query);
    ok($query->start,1);
    ok($query->end,35);
    ok($query->length,35);
    ok($query->dna,$alignments[0]->dna);
    ok($alignments[0]->strand,1);
    ok($query->strand,-1);

    ok($alignments[0]->get_tag_values('FLAGS'),$alignments[0]->flag_str);

    my @f = $sam->features(-name=>'EAS114_45:2:1:1140:1206');
    ok(scalar @f,2);

    @f   = $sam->features(-name=>'EAS114_45:2:1:1140:1206',
			  -attributes=>{FIRST_MATE=>1});
    ok(scalar @f,1);

    # try iteration
    my $i = $sam->get_seq_stream(
	-seq_id => 'seq2',
	-start  => 500,
	-end    => 800
	);
    ok($i);
    my $count = 0;
    while ($i->next_seq) { $count++ }
    ok($count,442);

    $i = $sam->get_seq_stream(); # all features!
    ok($i);
    $count = 0;
    while ($i->next_seq) { $count++ }
    ok($count,3307);

    # try the read_pair aggregation
    my @pairs = $sam->features(-type=>'read_pair',
			       -seq_id => 'seq2');
    ok(scalar @pairs,939);

    # try coverage
    my @coverage = $sam->features(-type   => 'coverage',
				  -seq_id => 'seq2');
    ok(scalar @coverage,1);
    my ($c) = $coverage[0]->get_tag_values('coverage');
    ok($c);
    ok($c->[0],3);
    ok($c->[1],4);
    ok($coverage[0]->type,"coverage:1584");

    exit 0;

# this is not a unit test, but a piece of demo to show cigar string
# processing
    for my $a (@alignments) {
	warn $a->display_name,' :: ',$a->flag_str,' :: ',
	$a->start,'..',$a->end,' ',' mpos=',$a->mate_start,' ins=',$a->isize,
	' qlen=',$a->cigar2qlen,
	' [',$a->strand,'/',$a->mstrand,']',
	' ',
	$a->cigar_str,
	' ',
	$a->mate_start,'..',$a->mate_end,
	"\n";
    }
}
