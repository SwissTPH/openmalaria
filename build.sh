#/bin/bash

# Assuming you have installed:
# gsl, git, cmake, xsd, xerces-c, boost

# Ubuntu:
# sudo apt-get install git cmake libboost-dev libgsl-dev libxerces-c-dev xsdcxx

# Mac OS:
# brew install git boost coreutils cmake gcc gsl xerces-c xsd

# Windows - Cygwin (MobaXTerm)
# apt-get install zlib-devel make python2 git cmake libboost-devel libgsl-devel xsd libxerces-c-devel

echo "=========================================================================="
echo "IMPORTANT: make sure dependencies are installed."
echo "Ubuntu/Debian:"
echo "sudo apt-get install build-essential git cmake libboost-dev libgsl-dev libxerces-c-dev xsdcxx"
echo "Mac OS:"
echo "brew install boost coreutils cmake gcc gsl xerces-c xsd"
echo "Windows - Cygwin (MobaXTerm):"
echo "apt-get install zlib-devel make python2 git cmake libboost-devel libgsl-devel xsd libxerces-c-devel"
echo "=========================================================================="

# stop on error
set -e

OMGIT=openmalaria
RELEASE=openMalaria
BRANCH=master

# For Windows - Cygwin (MobaXTerm)
export PATH=/usr/bin:$PATH

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE=Linux;;
    Darwin*)    MACHINE=Mac;;
    CYGWIN*)    MACHINE=Cygwin;;
    MINGW*)     MACHINE=MinGw;;
    *)          MACHINE="UNKNOWN:${unameOut}"
esac

echo "Running for ${MACHINE}"

# Clone git repo
if [ ! -d "$OMGIT" ] ; then
    git clone --branch $BRANCH https://github.com/SwissTPH/openmalaria.git $OMGIT
else
    pushd "$OMGIT" && git checkout $BRANCH && git pull && popd
fi

# Compile OpenMalaria
pushd $OMGIT/
mkdir -p build && pushd build && cmake -DCMAKE_BUILD_TYPE=Release -DOM_BOXTEST_ENABLE=OFF -DOM_CXXTEST_ENABLE=OFF .. && make -j4 && popd
popd

# Get version number
VERSION=$(cat $OMGIT/version.txt | cut -d'-' -f2)

# Prepare release folder by copying the necessary files
mkdir -p $RELEASE-$VERSION
cp $OMGIT/build/openMalaria $RELEASE-$VERSION/
cp -r $OMGIT/util/example/* $RELEASE-$VERSION/
cp $OMGIT/schema/scenario_41.xsd $RELEASE-$VERSION/
cp $OMGIT/test/densities.csv $RELEASE-$VERSION/
cp $OMGIT/test/autoRegressionParameters.csv $RELEASE-$VERSION/

# Compress into a tarball
tar -cavf $RELEASE-$VERSION.tar.gz $RELEASE-$VERSION/

# if Cygwin, copy dll files
if [ "${MACHINE}" = "Cygwin" ]; then
	pushd $RELEASE-$VERSION/
	mv $RELEASE-$VERSION/openMalaria $RELEASE-$VERSION/openMalaria.exe
	rm -f dlls; for i in $(ldd openMalaria.exe); do echo $i | grep "/usr" >> dlls;  done; for i in $(cat dlls); do cp -f $i .; done; rm -f dlls;
	popd
fi

echo "BUILD SUCCESS: $RELEASE-$VERSION.tar.gz"

# Don't delete the temp folders yet, let the user decide
# rm -rf $OMGIT $RELEASE-$VERSION/
