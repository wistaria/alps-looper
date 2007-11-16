#!/bin/sh

# for echo without trailing newline character
case `echo "testing\c"; echo 1,2,3`,`echo -n testing; echo 1,2,3` in
  *c*,-n*) ECHO_N= ECHO_C='
' ;;
  *c*,*  ) ECHO_N=-n ECHO_C= ;;
  *)       ECHO_N= ECHO_C='\c' ;;
esac

basedir=`dirname $0`
echo $ECHO_N `awk '$2=="LOOPER_DATE" {print $3}' $basedir/../looper/version.h | sed s/\"//g`$ECHO_C
