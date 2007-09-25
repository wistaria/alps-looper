#!/bin/sh

SRCDIR=`dirname $0`
PROG=$1; shift
BASE=$1; shift
PLOTS="$*"

if test "$SRCDIR" = "."; then
  SRCDIR="$PWD"
fi

if test -z "$NUM_PROCS"; then
  if test -f "/proc/cpuinfo"; then
    NUM_PROCS=`cat /proc/cpuinfo | grep processor | wc -l`
  else
    NUM_PROCS=1
  fi
fi

set -x

if test ! -f "$BASE.run/$BASE.out.xml"; then
  mkdir -p $BASE.run
  (cd $BASE.run && parameter2xml $SRCDIR/$BASE $BASE)
fi

if test `type lamnodes 1>&- 2>&-; echo $?` == 0; then
  if test `lamnodes > /dev/null 2>&1; echo $?` != 0; then
    lamboot;
  fi
fi

mpirun -x LD_LIBRARY_PATH -np $NUM_PROCS $PROG --evaluate $BASE.run/$BASE.in.xml

rm -f $BASE.db
archive --db-file=$BASE.db --command=install
archive --db-file=$BASE.db --command=append --xml-path=$BASE.run

for p in $PLOTS; do
  archive --db-file=$BASE.db --command=plot --plot-file=$SRCDIR/$p --output-path=.
  mv -f $p.text $BASE-$p.dat
done
