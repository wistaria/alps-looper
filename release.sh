#!/bin/sh

MAJOR=$1
MINOR=$2
REVISION=$3

if test -z "$REVISION"; then
    echo "$0 MAJOR MINOR REVISION"
    exit 127
fi

VERSION="$1.$2-r$REVISION"
echo "VERSION = $VERSION"

URL=`LANG=C svn info | grep 'URL:' | awk '{print $2}'`

DIR="alps-looper-$VERSION"
svn export -r $REVISION $URL $DIR

# set version
awk -v MAJOR=$MAJOR -v MINOR=$MINOR -v REV="r$REVISION" '$1=="set(LOOPER_VERSION_MAJOR" {$2=MAJOR")"} $1=="set(LOOPER_VERSION_MINOR" {$2=MINOR")"} $1=="set(LOOPER_VERSION_BUILD" {$2=REV")"} {print}' $DIR/CMakeLists.txt > $DIR/CMakeLists.txt.new
mv -f $DIR/CMakeLists.txt.new $DIR/CMakeLists.txt

# remove directories
(cd $DIR && rm -rf check extras/localsus extras/simple *.sh standalone v3.1)

# make tar.bz2
tar jcf $DIR.tar.bz2 $DIR
