#!/bin/sh

for p in aclocal automake autoconf
do
  ./missing $p
done
