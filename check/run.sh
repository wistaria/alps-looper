#!/bin/sh

SRCDIR=`dirname $0`
BASE=$1

case "$BASE" in
  *-d)
    REPL="diagonalization"
    ;;
  *-i)
    REPL="Ising"
    ;;
  *-p)
    REPL="path integral"
    ;;
  *-s)
    REPL="SSE"
    ;;
esac

if test ! -f "$BASE.run/$BASE.out.xml"; then
    echo "REPRESENTATION = \"$REPL\"" > $BASE
    cat $SRCDIR/site >> $BASE
    mkdir -p $BASE.run
    (cd $BASE.run && parameter2xml ../$BASE $BASE && cp -f $BASE.in.xml $BASE.out.xml)
fi
../loop --Tmin 1 $BASE.run/$BASE.out.xml
../loop_evaluate $BASE.run/$BASE.task*.out.xml
rm -f $BASE.db
archive --db-file=$BASE.db --command=install
for x in `ls $BASE.run/$BASE.task*.out.xml`; do
    archive --db-file=$BASE.db --command=append --xml-path=$x
done
