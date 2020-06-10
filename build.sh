#/bin/bash

# Assuming you have installed:
# gsl, git, cmake, xsd, xerces-c, boost

# stop on error
set -e

# Clone to
OMGIT=openmalaria

# Clone from
RELEASE=openMalaria
BRANCH=master

# option
CLEAN=0
TESTS=OFF
JOBS=4
ARTIFACT=

# shell
CYGWIN=0

# For Windows - Cygwin (MobaXTerm)
export PATH=/usr/bin:$PATH

function printHelp {
    echo "HELP: make sure dependencies are installed (adapt to your distribution)"
    echo "======================================================================="
    echo "Ubuntu/Debian: sudo apt-get install build-essential git cmake libboost-dev libgsl-dev libxerces-c-dev xsdcxx"
    echo "--------------"
    echo "Mac OS: brew install git boost coreutils cmake gcc gsl xerces-c xsd"
    echo "--------------"
    echo "Windows - Cygwin (MobaXTerm): apt-get install gcc-g++ git cmake make python2 zlib-devel libboost-devel libgsl-devel xsd libxerces-c-devel"
}

function parseArguments {
    for i in "$@"; do
        case $i in
            -c|--clean)         CLEAN=1 && shift ;;
            -t|--tests)         TESTS=ON && shift ;;
            -a=*|--artifact=*)  ARTIFACT=${i#*=} && shift ;;
            -j=*|--jobs=*)      JOBS="${i#*=}" && shift ;;
            -h|--help)          printHelp && shift && exit ;;
            *) ;;
        esac
    done
}

function isWindows {
    unameOut="$(uname -s)"
    case "${unameOut}" in
        CYGWIN*)    CYGWIN=1;;
        MINGW*)     CYGWIN=1;;
        *)          CYGWIN=0;;
    esac
}

function clone {
    # Already in openmalaria?
    if [ -d .git ]; then
        if [ $(basename -s .git `git remote get-url origin`) = "openmalaria" ]; then
            echo "Already in openmalaria repo, not cloning."
        else
            echo "Error: this git repository is not openmalaria."
        fi
    else
        # Clone git repo
        if [ ! -d "$OMGIT" ] ; then
            git clone --branch $BRANCH https://github.com/SwissTPH/openmalaria.git $OMGIT
        else
            echo "Folder $OMGIT already exist, not cloning."
            cd "$OMGIT" && git checkout $BRANCH && git pull
        fi
    fi
}

function build {
    mkdir -p build

    if [ $CLEAN -eq 1 ]; then
        pushd build && rm -rf * && popd
    fi

    # Compile OpenMalaria
    pushd build
    cmake -DCMAKE_BUILD_TYPE=Release -DOM_BOXTEST_ENABLE=$TESTS -DOM_CXXTEST_ENABLE=$TESTS .. && make -j$JOBS
    popd
}

function test {
    if [ $TESTS = "ON" ]; then
        pushd build && ctest -j$JOBS && popd
    fi
}

function package {
    # Get version number
    VERSION=$(cat version.txt | cut -d'-' -f2)
    

    echo "ARTIFACT IS $ARTIFACT"

    if [ -z "${ARTIFACT}" ]; then
        ARTIFACT=$RELEASE-$VERSION
        echo $ARTIFACT
    fi

    # Prepare release folder by copying the necessary files
    mkdir -p $ARTIFACT
    cp build/openMalaria $ARTIFACT/
    cp -r util/example/* $ARTIFACT/
    cp schema/scenario_41.xsd $ARTIFACT/
    cp test/densities.csv $ARTIFACT/
    cp test/autoRegressionParameters.csv $ARTIFACT/

    # if Cygwin, copy dll files
    if [ $CYGWIN -eq 1 ]; then
        pushd $ARTIFACT/
        #mv openMalaria openMalaria.exe
        rm -f dlls
        ldd openMalaria
        for i in $(ldd openMalaria); do
            echo $i
            echo $i | grep "/usr" >> dlls
        done
        cat dlls
        for i in $(cat dlls); do
            cp -f $i .
        done
        rm -f dlls
        popd
    fi

    # Compress
    if [ $CYGWIN -eq 1 ]; then
        7z a $ARTIFACT.zip "$ARTIFACT/"
    else
        tar -cavf $ARTIFACT.tar.gz $ARTIFACT/
    fi
}

parseArguments $@ # CLEAN=0:1, JOBS=1:N, TESTS=ON:OFF, print help
isWindows # set CYGWIN=0:1
clone # clone the repo
build # compile
test # test
package # create $ARTIFACT.zip or $ARTIFACT.tar.gz

echo "BUILD SUCCESS: $(ls $ARTIFACT.*)"

# Don't delete the temp folders yet, let the user decide
# rm -rf $OMGIT $ARTIFACT/
