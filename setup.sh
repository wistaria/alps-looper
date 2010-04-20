#!/bin/sh

set -x

DIR="$HOME/development/alps"
if test -d "$DIR"; then
  cp -rp $DIR/config/update_preamble config/
fi
