#!/bin/sh

VERSION=$1
if test -z "$VERSION"; then
    echo "$0 VERSION"
    exit 127
fi
echo "VERSION = $VERSION"
DIR="alps-looper-$VERSION"

set -x

git archive --prefix=$DIR/ --format=tar $VERSION | bzip2 > $DIR.tar.bz2
git archive --prefix=$DIR/ --format=zip $VERSION > $DIR.zip
