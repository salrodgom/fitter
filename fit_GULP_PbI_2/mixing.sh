#!/bin/bash
while read line ; do
 atom_types=$(echo $line | awk '{print $1,$2}')
 eii=$(echo $line | awk '{print $3}')
 rii=$(echo $line | awk '{print $4}')
 sigma=$(echo "scale=4; 0.890898718*${rii}" | bc -lq)
 scale=$(echo $line | awk '{print $5}')
 bineii=$(./real2binary_string_IEEE $eii)
 binrii=$(./real2binary_string_IEEE $sigma)
 echo "$bineii # $eii"
 echo "$binrii # $sigma"
done < mixing_uff.txt
