#!/bin/sh

dir=`dirname $0`

touch $dir/aclocal.m4
find $dir -name Makefile.in | xargs touch
touch $dir/configure
find $dir -name config.h.in | xargs touch
