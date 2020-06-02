#/bin/bash

# Assuming you have installed:
# gsl, git, xsd, xerces-c, boost

# Ubuntu:
# sudo apt-get install git libboost-dev libgsl-dev libxerces-c-dev xsdcxx  

echo "IMPORTANT: make sure you have all the dependencies installed."
echo "On Ubuntu you can install them with:"
echo "sudo apt-get install git libboost-dev libgsl-dev libxerces-c-dev xsdcxx"

# stop on error
set -e

OMGIT=openmalaria
RELEASE=openMalaria

git clone https://github.com/SwissTPH/openmalaria.git $OMGIT

pushd $OMGIT/
mkdir -p build && pushd build && cmake -DCMAKE_BUILD_TYPE=Release -DOM_BOXTEST_ENABLE=OFF -DOM_CXXTEST_ENABLE=OFF .. && make -j4 && popd
VERSION=$(cat version.txt | cut -d'-' -f2)
popd

mkdir -p $RELEASE-$VERSION
cp $OMGIT/build/openMalaria $RELEASE-$VERSION/
cp -r $OMGIT/util/example/* $RELEASE-$VERSION/
cp $OMGIT/schema/scenario_41.xsd $RELEASE-$VERSION/
cp $OMGIT/test/densities.csv $RELEASE-$VERSION/
cp $OMGIT/test/autoRegressionParameters.csv $RELEASE-$VERSION/

tar -cavf $RELEASE-$VERSION.tar.gz $RELEASE-$VERSION/

echo "BUILD SUCCESS: $RELEASE-$VERSION.tar.gz"

#echo "You can remove $OMGIT and $RELEASE-$VERSION folders"
# deleting temporary working files
# rm -rf $OMGIT $RELEASE-$VERSION/
