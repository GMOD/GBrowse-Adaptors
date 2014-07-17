use Module::Build;
# See perldoc Module::Build for details of how this works

Module::Build->new
    ( module_name     => 'Bio::DB::Das::BioSQL',
      license         => 'perl',
      requires        => {
                           'Bio::DB::BioDB'     => 1.006,
                           'Bio::Root::Version' => 1.006,
                         }
    )->create_build_script;
