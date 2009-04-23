# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 6;

BEGIN { 
  use_ok( 'Bio::DB::Das::BioSQL' ); 
  use_ok( 'Bio::DB::Das::BioSQL::BioDatabaseAdaptor' ); 
  use_ok( 'Bio::DB::Das::BioSQL::DBAdaptor' ); 
  use_ok( 'Bio::DB::Das::BioSQL::Iterator' ); 
  use_ok( 'Bio::DB::Das::BioSQL::PartialSeqAdaptor' ); 
  use_ok( 'Bio::DB::Das::BioSQL::Segment' );
}



