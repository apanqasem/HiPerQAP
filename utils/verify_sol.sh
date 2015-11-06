#!/bin/sh

prog=$1
datafile=$2
sols=$3
seed=$4

solfile=sol.txt

[ `which costcheck` ] || { echo "could not find costcheck. make sure QAP bin directory is in your path"; exit 0; }

[ -x src ] || { echo "verify_sol.sh must be run in parent of src directory"; exit 0;}


cd src
./$prog $datafile $sols $seed | grep -E "size|best solution|best cost" | awk -F ":" '{print $2}' > $solfile

costcheck $datafile $solfile 

rm -rf $solfile

cd ../


