#!/bin/sh

# find . -name Makefile.in | xargs ./replace_header.pl
# find . -name '*.mk.in' | xargs ./replace_header.pl
# find . -name '*.mk' | xargs ./replace_header.pl
find . -name '*.h.in' | xargs ./replace_header.pl
find . -name '*.h' | xargs ./replace_header.pl
find . -name '*.C.in' | xargs ./replace_header.pl
find . -name '*.C' | xargs ./replace_header.pl

find . -name '*.orig' | xargs rm
