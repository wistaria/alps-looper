#!/bin/sh

dir=`dirname $0`

touch $dir/aclocal.m4
find $dir -name Makefile.in | xargs touch
touch $dir/configure
touch $dir/looper/config.h.in
