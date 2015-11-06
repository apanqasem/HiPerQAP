#!/bin/bash

datasets="tai30a tai30b tai40a tai40b tai50a tai50b tai60a tai60b tai80a tai80b tai100a tai100b lipa70a lipa90"
progs="tabu_dyn"

for prog in $progs; do
  for data in $datasets; do
      qap.sh $prog $data 1024 1 1
  done
done

