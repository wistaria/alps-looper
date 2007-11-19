#!/bin/sh

# for echo without trailing newline character
case `echo "testing\c"; echo 1,2,3`,`echo -n testing; echo 1,2,3` in
  *c*,-n*) ECHO_N= ECHO_C='
' ;;
  *c*,*  ) ECHO_N=-n ECHO_C= ;;
  *)       ECHO_N= ECHO_C='\c' ;;
esac

VERSION_H=`find . -name version.h`

echo $ECHO_N `awk 'NF==3 && $2 ~ /_DATE$/ {print $3}' $VERSION_H | sed s/\"//g`$ECHO_C
