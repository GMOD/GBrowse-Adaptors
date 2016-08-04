#!/bin/bash

set -uex

# for Bio::DB::BigWig
SOURCE_KENTSRC="ftp://ftp.sanger.ac.uk/pub/cancer/legacy-dependancies/jksrc.v334.zip"
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

#"File::ShareDir" "File::ShareDir::Install" "Const::Fast" "File::Which"
perlmods=( "ExtUtils::CBuilder" "Module::Build~0.42" "LWP::UserAgent" "Bio::Root::Version~1.006009001")
for i in "${perlmods[@]}" ; do
  $CPANM -v --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

cd $SETUP_DIR
rm -rf kent kentsrc.zip
get_distro "kentsrc" $SOURCE_KENTSRC
unzip -q kentsrc.zip
perl -pi -e 's/(\s+CFLAGS=)$/${1}-fPIC/' kent/src/inc/common.mk
cd kent/src/lib
export MACHTYPE=i686    # for a 64-bit system
make -j$CPU
cd ../
export KENT_SRC=`pwd`
cd $SETUP_DIR
mkdir bigfile
get_distro "bigfile" $SOURCE_BIGFILE
tar --strip-components 1 -C bigfile -zxf bigfile.tar.gz
cd bigfile
perl Build.PL --install_base=$INST_PATH
./Build
./Build test
./Build install
rm -f kentsrc.zip


