# *Nix shell script compiles creates debug and release build directories
# Takes the following optional arguments:
# -debugdir: Path to debug build directory
# -releasedir: Path to release build directory
# -builddir: Path to custom build directory (something other than debug|release)
# <target>: The name of the build target 

# Initialize to default paths
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
LOOPS_DEBUG_DIR=$SCRIPT_DIR/../CMakeBuild/Debug
LOOPS_RELEASE_DIR=$SCRIPT_DIR/../CMakeBuild/Release
# If build target is not 'all', 'debug', or 'release'
LOOPS_BUILD_DIR=
LOOPS_BUILD_TARGET=all

# Function to run make in parameter build directory
function doMake {
    LOOPS_MAKE_BUILD_DIR=$1
    echo Running make for target dir $LOOPS_MAKE_BUILD_DIR
    # Navigate to target build directory
    pushd .
    cd $LOOPS_MAKE_BUILD_DIR
    # Run make in target build directory
    make
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
    -builddir)
        LOOPS_BUILD_DIR=$2
        shift
        ;;
    *)
        LOOPS_BUILD_TARGET=$1
        ;;
    esac
    shift
done

# Determine build directories
if [ -n "$LOOPS_BUILD_DIR" ]
then
    echo Building custom target at $LOOPS_BUILD_DIR
    doMake $LOOPS_BUILD_DIR
elif [ "$LOOPS_BUILD_TARGET" == "all" ]
then
    echo Building all targets
    doMake $LOOPS_DEBUG_DIR
    doMake $LOOPS_RELEASE_DIR
elif [ "$LOOPS_BUILD_TARGET" == "debug" ]
then
    echo Building debug target at $LOOPS_DEBUG_DIR
    doMake $LOOPS_DEBUG_DIR
elif [ "$LOOPS_BUILD_TARGET" == "release" ]
then
    echo Building release target at $LOOPS_RELEASE_DIR
    doMake $LOOPS_RELEASE_DIR
else
    echo No build found
fi

