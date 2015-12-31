#!/bin/bash

data=$1
cost=$2

tai30a=1818146
tai30b=637117113
tai35a=2422002
tai35b=283315445
tai40a=3139370
tai40b=637250948
tai50a=4938796
tai50b=458821517
tai60a=7205962
tai60b=608215054
tai80a=13515450
tai80b=818415043
tai100a=21054656
tai100b=1185996137
lipa70a=169755
lipa90=360630


[ $data = "tai35a" ] && { echo $tai35a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'; exit 0; }
[ $data = "tai35b" ] && { echo $tai35b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'; exit 0; }

if [ $data = "tai30a" ]; then 
   echo $tai30a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai30b" ]; then 
   echo $tai30b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai40a" ]; then 
   echo $tai40a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai40b" ]; then 
   echo $tai40b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai50a" ]; then 
   echo $tai50a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai50b" ]; then 
   echo $tai50b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai60a" ]; then 
   echo $tai60a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai60b" ]; then 
   echo $tai60b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai80a" ]; then 
   echo $tai80a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai80b" ]; then 
   echo $tai80b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai100a" ]; then 
   echo $tai100a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "tai100b" ]; then 
   echo $tai100b $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "lipa70a" ]; then 
   echo $lipa70a $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
if [ $data = "lipa90" ]; then 
   echo $lipa90 $cost | awk '{ printf "%3.2f%%\n", (($2 - $1)/$1 * 100) }'
fi
