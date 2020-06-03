#/bin/bash

# Assuming you have installed:
# gsl, git, xsd, xerces-c, boost

# Ubuntu:
# sudo apt-get install git libboost-dev libgsl-dev libxerces-c-dev xsdcxx

brew install boost@1.60 coreutils cmake gcc gsl xerces-c xsd

echo "=========================================================================="
echo "IMPORTANT: make sure dependencies are installed."
echo "Ubuntu/Debian:"
echo "sudo apt-get install build-essential git cmake libboost-dev libgsl-dev libxerces-c-dev xsdcxx"
echo "Mac OS:"
echo "brew install boost coreutils cmake gcc gsl xerces-c xsd"
echo "=========================================================================="

# stop on error
set -e

OMGIT=openmalaria
RELEASE=openMalaria
BRANCH=master

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

echo "BUILD SUCCESS: $RELEASE-$VERSION.tar.gz"

# Don't delete the temp folders yet, let the user decide
# rm -rf $OMGIT $RELEASE-$VERSION/
