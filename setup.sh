#!/bin/sh

set -x

DIR="$HOME/development/alps"
if test -d "$DIR"; then
  cp -rp $DIR/bootstrap .
  cp -rp $DIR/config/ac_alps.m4 config/
  cp -rp $DIR/config/run-test config/
  cp -rp $DIR/config/update_preamble config/
fi

DIR="$HOME/development/simtool"
if test -d "$DIR"; then
  cp -fp $DIR/config/extract_date.sh config/
  cp -fp $DIR/config/extract_version.sh config/
  cp -fp $DIR/config/update_version config/
fi
