#!/bin/bash 
declare -A e
declare -A s
declare -A atom_types
i=0
while read line ; do
 let i++
 #bineii=$(./real2binary_string_IEEE $eii)
 #binrii=$(./real2binary_string_IEEE $sigma)
 #echo "$bineii # $eii   $atom_types"
 #echo "$binrii # $sigma $atom_types"
 e[${i},${i}]=$(echo $line | awk '{print $3}')
 rii=$(echo $line | awk '{print $4}')
 s[${i},${i}]=$(echo "scale=4; 0.890898718*${rii}" | bc -lq)
 scale=$(echo $line | awk '{print $5}')
 atom_types[${i}]=$(echo $line | awk '{print $1}')
done < mixing_uff.txt
n_atoms=$i
for ((i=1;i<=n_atoms;i++)) do
 for ((j=i;j<=n_atoms;j++)) do
  e[${i},${j}]=$(echo "scale=6; e(0.5*l(${e[$i,$i]} * ${e[${j},${j}]}))" | bc -lq)
  s[$i,$j]=$(echo "scale=6; 0.5*(${s[$i,$i]}+${s[$j,$j]})" | bc -lq)
  echo ${e[$i,$j]} ${s[$i,$j]} ${atom_types[$i]} ${atom_types[$j]}
 done
done
