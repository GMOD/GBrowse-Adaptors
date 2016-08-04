#!/bin/bash

set -uex

# for Bio::DB::BigWig
#SOURCE_KENTSRC="ftp://ftp.sanger.ac.uk/pub/cancer/legacy-dependancies/jksrc.v334.zip"
SOURCE_KENTSRC="ftp://ftp.sanger.ac.uk/pub/cancer/legacy-dependancies/jksrc.v336.zip"
SOURCE_BIGFILE="http://www.cpan.org/authors/id/L/LD/LDS/Bio-BigFile-1.07.tar.gz"

CPU=`grep -c ^processor /proc/cpuinfo`

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

INST_PATH=$1
if [ "$#" -eq "2" ] ; then
  SOURCE_KENTSRC=$2
fi

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH
cd $INST_PATH
INST_PATH=`pwd`
mkdir -p $INST_PATH/bin
cd $INIT_DIR

unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
export PATH="$INST_PATH/bin:$PATH"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`

if [ -e $SETUP_DIR/basePerlDeps.success ]; then
  echo "Previously installed base perl deps..."
else
  perlmods=( "ExtUtils::CBuilder" "Module::Build~0.42" "LWP::UserAgent" "Bio::Root::Version~1.006009001")
  for i in "${perlmods[@]}" ; do
    $CPANM -v --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
  done
  touch $SETUP_DIR/basePerlDeps.success
fi

cd $SETUP_DIR
if [ -e $SETUP_DIR/kentsrc.success ]; then
  echo "Previously unpacked kentsrc..."
else
  rm -rf kent kentsrc.zip
  get_distro "kentsrc" $SOURCE_KENTSRC
  unzip -q kentsrc.zip
  perl -pi -e 's/(\s+CFLAGS=)$/${1}-fPIC/' kent/src/inc/common.mk
  perl -pi -e 'if($_ =~ m/^CFLAGS/ && $_ !~ m/\-fPIC/i){chomp; s/#.+//; $_ .= " -fPIC -Wno-unused -Wno-unused-result\n"};' kent/src/htslib/Makefile
  touch $SETUP_DIR/kentsrc.success
fi

cd kent/src/htslib
make -j$CPU
cd ../lib
export MACHTYPE=i686    # for a 64-bit system
make -j$CPU
cd ../
export KENT_SRC=`pwd`
cd $SETUP_DIR
cd $INIT_DIR/Bio-BigFile
perl Build.PL --install_base=$INST_PATH
./Build clean
./Build
./Build test
./Build install
rm -f kentsrc.zip


