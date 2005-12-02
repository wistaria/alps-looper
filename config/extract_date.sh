#! /bin/sh

basedir=`dirname $0`
echo -n `awk '$2=="LOOPER_DATE" {print $3}' $basedir/../looper/version.h | sed s/\"//g`
