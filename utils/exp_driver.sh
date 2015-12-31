#!/bin/bash

datasets="tai30a tai30b tai40a tai40b tai50a tai50b tai60a tai60b tai80a tai80b tai100a tai100b lipa70a lipa90"
datasets="tai40a tai40b"
progs="tabu_dyn"

# for prog in $progs; do
#   for data in $datasets; do
#       qap.sh $prog $data 1024 40 40
#   done
# done

instances="32 64 96 128 256 512 768 1024"
restarts="1 5 10 15 20 25 30 35 40"

for prog in $progs; do
	for data in $datasets; do
#		for hoods in {1..40}; do
#		for inst in ${instances}; do
		for r in ${restarts}; do
      qap.sh $prog $data 256 10 $r
		done
#	done
		done
done

