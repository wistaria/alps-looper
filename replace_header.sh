#!/bin/sh

find . -name Makefile.in | xargs ./replace_header.pl
find . -name '*.h' | xargs ./replace_header.pl
find . -name '*.C' | xargs ./replace_header.pl
