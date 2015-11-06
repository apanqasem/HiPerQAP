#!/bin/bash

#
# Executes an instance of qap with supplied number of initial solutions, datafile, number of neighborhoods
# and number of random restarts
#

if [ $# -ne 5 ]; then                                                          
  echo "usage :"                                                               
  echo "    run_exp_single.sh prog datafile solns hoods restarts"
  echo "       prog = executable"
  echo "       solns = number of solutions"
  echo "       datafile = name of dataset (path and extension added automatically)"
  echo "       hoods = number of neighborhoods to explore"
  echo "       restarts = number of random restarts"
  exit                                                                         
fi                                                                             

prog=$1
datafile=$2
solns=$3
hoods=$4
restarts=$5

outpath=`pwd`
DATAPATH=../datasets
outfile=${outpath}/perf.dat


[ "$QAPHOME" ] || { echo "QAPHOME env var not defined. Set QAPHOME to install dir"; exit 0;}

rm -rf ${outfile}

srcdir=${QAPHOME}/src

[ -x ${srcdir} ] || { echo "could not find src directory"; exit 0;}

pushd ${srcdir} &> /dev/null
i=0
while [ $i -lt ${restarts} ]; do 
  RANDSEED=$RANDOM 
  echo -n $prog $datafile $RANDSEED " " >> ${outfile}
  vals=`./$prog $DATAPATH/$datafile.txt $solns $RANDSEED $hoods\
        | grep -E "cost|kernel" | awk '{printf $3 " "}'`

  cost=`echo $vals | awk '{print $1}'`
  time=`echo $vals | awk '{print $2}'`
  gap=`../utils/calc_gap.sh $datafile $cost`

  echo $gap $time >> ${outfile}

	i=$(($i+1))
done

popd &> /dev/null

