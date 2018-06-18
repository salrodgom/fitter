#!/bin/bash
# 0.0
# 0.179046289062E+05
n_refits=0
for x1 in 0.25 0.50 0.90 1.0 1.10 1.50 ; do              # 5
 x01=0.482217
 xx1=$(echo "$x01*$x1" | bc -lq)
 #
 for x2 in 0.50 0.90 1.0 1.50 2.0 ; do             # 5
  x02=17904.6289062
  xx2=$(echo "$x02*$x2" | bc -lq)
  for x3 in 2.0 20.0 200.0 700.0 ; do             # 5
   let n_refits++
   x03=0.0
   xx3=$(echo "$x03 + ($x3)" | bc -lq)
   #
   xxx1=$(./real2binary_string_IEEE $xx1)
   xxx2=$(./real2binary_string_IEEE $xx2)
   xxx3=$(./real2binary_string_IEEE $xx3)
   echo $xxx1 "#" $xx1 $n_refits
   echo $xxx2 "#" $xx2
   echo $xxx3 "#" $xx3
  done
 done
done
