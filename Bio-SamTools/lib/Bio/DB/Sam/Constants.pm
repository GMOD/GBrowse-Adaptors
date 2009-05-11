use Bio::DB::Sam::Constants;

use strict;

BEGIN {
    require Exporter;
    our @ISA    = qw(Exporter);
    our @EXPORT = qw(CIGAR_SYMBOLS BAM_CIGAR_SHIFT BAM_CIGAR_MASK
                 BAM_CMATCH BAM_CINS BAM_CDEL BAM_CREF_SKIP
                 BAM_CSOFT_CLIP BAM_CHARD_CLIP BAM_CPAD FLAGS);
    our @EXPORT_OK = (@EXPORT,'RFLAGS');
}

use constant CIGAR_SYMBOLS   => [qw(M I D N S H P)];
use constant BAM_CIGAR_SHIFT => 4;
use constant BAM_CIGAR_MASK  => (1 << BAM_CIGAR_SHIFT) - 1;
use constant BAM_CMATCH      => 0;
use constant BAM_CINS        => 1;
use constant BAM_CDEL        => 2;
use constant BAM_CREF_SKIP   => 3;
use constant BAM_CSOFT_CLIP  => 4;
use constant BAM_CHARD_CLIP  => 5;
use constant BAM_CPAD        => 6;

use constant FLAGS => {
    0x0001 => 'PAIRED',
    0x0002 => 'MAP_PAIR',
    0x0004 => 'UNMAPPED',
    0x0008 => 'M_UNMAPPED',
    0x0010 => 'REVERSED',
    0x0020 => 'M_REVERSED',
    0x0040 => 'FIRST_MATE',
    0x0080 => 'SECOND_MATE',
    0x0100 => 'NOT_PRIMARY',
    0x0200 => 'QC_FAILED',
    0x0400 => 'DUPLICATE'
};
use constant RFLAGS => {reverse %{FLAGS()}};

1;
