#!/bin/bash
# 0.0
# 0.179046289062E+05
n_refits=0
for x1 in 0.6 1.0 ; do              # 5
 x01=70359906.629702
 xx1=$(echo "$x01*$x1" | bc -lq)
 for x2 in 1.0 ; do             # 5
  x02=0.131258
  xx2=$(echo "$x02*$x2" | bc -lq)
  for x3 in 0.0 ; do             # 5
   x03=0.0
   xx3=$(echo "$x03 + ($x3)" | bc -lq)
   for x4 in 0.60 1.0 1.0 ; do              # 5
    x04=103496.133010
    xx4=$(echo "$x04*$x4" | bc -lq)
    for x5 in 0.60 1.0 1.2 ; do              # 5
     x05=0.321737
     xx5=$(echo "$x05*$x5" | bc -lq)
     for x6 in 0.0 ; do              # 5
      x06=0.0
      xx6=$(echo "$x06*$x6" | bc -lq)
      for x7 in 0.9 1.0 1.2 ; do              # 5
       x07=22793.338582
       xx7=$(echo "$x07*$x7" | bc -lq)
       for x8 in 0.60 1.0 1.2 ; do              # 5
        x08=0.482217
        xx8=$(echo "$x08*$x8" | bc -lq)
        for x9 in 0.0 ; do              # 5
         x09=696.949542
         xx9=$(echo "$x09*$x9" | bc -lq)
         #
         let n_refits++
         xxx1=$(./real2binary_string_IEEE $xx1)
         xxx2=$(./real2binary_string_IEEE $xx2)
         xxx3=$(./real2binary_string_IEEE $xx3)
         xxx4=$(./real2binary_string_IEEE $xx4)
         xxx5=$(./real2binary_string_IEEE $xx5)
         xxx6=$(./real2binary_string_IEEE $xx6)
         xxx7=$(./real2binary_string_IEEE $xx7)
         xxx8=$(./real2binary_string_IEEE $xx8)
         xxx9=$(./real2binary_string_IEEE $xx9)
         echo $xxx1 "#" $xx1 $n_refits
         echo $xxx2 "#" $xx2
         echo $xxx3 "#" $xx3
         echo $xxx4 "#" $xx4
         echo $xxx5 "#" $xx5
         echo $xxx6 "#" $xx6
         echo $xxx7 "#" $xx7
         echo $xxx8 "#" $xx8
         echo $xxx9 "#" $xx9
         #
        done
       done
      done
     done
    done
   done
  done
 done
done
