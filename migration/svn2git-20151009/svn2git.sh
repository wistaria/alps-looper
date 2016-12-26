#!/bin/sh

REVFILE=alps-looper-rev.dat
DIR=$HOME/alps-looper
COMMITLOG=/tmp/commit.$$

REVS=$(cat $REVFILE)

for rev in $REVS; do
    echo [processing $rev]
    (cd $DIR; svn log -r$rev > $COMMITLOG)
    USER=$(awk '$2=="|" && $4=="|" && $12=="|" {print $3}' $COMMITLOG)
    DATE=$(awk '$2=="|" && $4=="|" && $12=="|" {print $5,$6,$7}' $COMMITLOG)
    echo "USER = $USER, DATE=$DATE"
    (cd $DIR; svn update -r$rev)
    REMOVED=$(cd $DIR; git status | awk '$1=="deleted:" {print $2}')
    for r in $REMOVED; do
	(cd $DIR; git rm $r)
    done
    (cd $DIR; git add *)
    (cd $DIR; git commit --file $COMMITLOG --date "$DATE")
done
