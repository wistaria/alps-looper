#!/bin/sh

BINDIR=$1

do_eval()
{
for p in $parameters
do
  echo "parameter = $p"

  # echo -n "Compacting "
  # for f in `ls $p.task[0-9]*.out.run[0-9]*`; do
  #   echo -n "."
  #   $BINDIR/compactrun $f > /dev/null
  # done
  # echo

  xmls=`ls $p.task[0-9]*.out.xml 2> /dev/null`

  if test -n "$xmls"; then
    echo -n "Evaluating "
    for f in $xmls; do
      echo -n "."
      ./../loop_evaluate $f > /dev/null 2>&1
    done
    echo

    $BINDIR/archivecat $xmls \
      | sed 's/Energy Density/ene/g' \
      | sed 's/Specific Heat/sh/g' \
      | sed 's/Staggered Magnetization\^2/smag2/g' \
      | sed 's/Staggered Magnetization\^4/smag4/g' \
      | sed 's/Staggered Magnetization/smag/g' \
      | sed 's/Staggered Susceptibility/ssus/g' \
      | sed 's/Magnetization\^2/umag2/g' \
      | sed 's/Magnetization\^4/umag4/g' \
      | sed 's/Magnetization/umag/g' \
      | sed 's/Susceptibility/usus/g' \
      > archive.xml
  
    file="$p.dat"
    rm -f $file

    for m in ene sh umag umag2 umag4 usus smag smag2 smag4 ssus; do
    echo "measurement = $m"

      cat <<EOF > plot.xml
<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="http://xml.comp-phys.org/2003/4/plot2html.xsl"?>
<plot>
  <xaxis label="T" type="PARAMETER" name="T"/>
  <yaxis label="$m" type="SCALAR_AVERAGE"/>  
</plot>
EOF
    $BINDIR/extracttext plot.xml archive.xml > plot-1.dat 2> /dev/null
    if test `wc plot-1.dat | awk '{print $1}'` == 0; then
      rm -f plot-1.dat
    else
      if test -f "$file"; then
        mv $file plot-0.dat
        join plot-0.dat plot-1.dat > $file
        rm -f plot-0.dat plot-1.dat
      else
        mv plot-1.dat $file
      fi
    fi

    rm -f plot.xml

    done
  fi
done
return 0
}

do_eval_g()
{
  for p in $parameters; do
  if test -f "$p.out"; then
    awk '$1=="T" {temp=$3} $1=="energy" && $2=="density" {ene=$4} $1=="specific" {sh=$4} $1=="uniform" && $2=="magnetization" && $3=="=" {umag=$4} $1=="uniform" && $2=="magnetization^2" {umag2=$4} $1=="uniform" && $2=="magnetization^4" {umag4=$4} $1=="uniform" && $2=="susceptibility" {usus=$4} $1=="staggered" && $2=="magnetization" && $3=="=" {smag=$4} $1=="staggered" && $2=="magnetization^2" {smag2=$4} $1=="staggered" && $2=="magnetization^4" {smag4=$4} $1=="staggered" && $2=="susceptibility" {print temp,ene,sh,umag,umag2,umag4,usus,smag,smag2,smag4,$4}' $p.out > $p.dat
  fi
  done
}

parameters='site-p site-e'
do_eval

parameters='site-g'
do_eval_g
