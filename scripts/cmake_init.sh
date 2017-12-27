#!/bin/bash 
# *Nix shell script which creates debug and release build
# directories using CMake
# Takes the following optional arguments:
# -debugdir: Path to debug build directory
# -releasedir: Path to release build directory
# -useintel: toggle usage of intel compiler

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
LOOPS_CODE_DIR=$SCRIPT_DIR/..
LOOPS_DEBUG_DIR=$LOOPS_CODE_DIR/CMakeBuild/Debug
LOOPS_RELEASE_DIR=$LOOPS_CODE_DIR/CMakeBuild/Release
USE_INTEL=0

# Function to perform an out-of-source cmake build
function doOutOfSourceBuild {
    LOOPS_BUILD_DIR=$1
    LOOPS_BUILD_TYPE=$2
    USE_INTEL_COMPILER=$3
    echo Generating out-of-source $LOOPS_BUILD_TYPE build at $LOOPS_BUILD_DIR using $LOOPS_CODE_DIR
    # Clear out any previous build directory
    rm -r $LOOPS_BUILD_DIR
    # Navigate to build directory
    mkdir -p $LOOPS_BUILD_DIR
    pushd .
    cd $LOOPS_BUILD_DIR
    if [ $USE_INTEL_COMPILER -eq 1 ]; then
        echo Using intel compiler
        cmake -DCMAKE_BUILD_TYPE=$LOOPS_BUILD_TYPE -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -G"CodeBlocks - Unix Makefiles" $LOOPS_CODE_DIR
    else
        echo Using default compiler
        cmake -DCMAKE_BUILD_TYPE=$LOOPS_BUILD_TYPE -G"CodeBlocks - Unix Makefiles" $LOOPS_CODE_DIR
    fi
    popd
}

# Check if there are any parameter overrides
while [ $# -gt 0 ]
do
    case "$1" in
    -debugdir)
        LOOPS_DEBUG_DIR=$2
        shift
        ;;
    -releasedir)
        LOOPS_RELEASE_DIR=$2
        shift
        ;;
    -useintel)
        USE_INTEL=1
        ;;
    esac
    shift
done

# Create debug build
doOutOfSourceBuild $LOOPS_DEBUG_DIR Debug $USE_INTEL

# Create release build
doOutOfSourceBuild $LOOPS_RELEASE_DIR Release $USE_INTEL
