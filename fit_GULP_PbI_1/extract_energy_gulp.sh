#!/bin/bash
if [ -f cell_energy_expansion_gulp.txt ] ; then rm cell_energy_expansion_gulp.txt ; touch cell_energy_expansion_gulp.txt ; fi
for scale in 1.1000 1.0500 0.9750 0.9625 1.0375 1.0125 1.0000 0.9875 0.9500 0.9000 0.9125 0.9250 1.0625 1.0750 0.8875 1.1125 ; do
 structure=MAPbI-${scale} 
 folder=${structure}_dir
 if [ -d $folder ] ; then
  energy=$(grep "Total lattice energy       =" ${folder}/${structure}.gout | grep "eV" | awk '{print $5}' | tail -n1 )
  gnorm=0.0 
  a=$(grep "_cell_length_a" ${folder}/test.cif | awk '{print $2}')
  b=$(grep "_cell_length_b" ${folder}/test.cif | awk '{print $2}')
  c=$(grep "_cell_length_c" ${folder}/test.cif | awk '{print $2}')
  alpha=$(grep "_cell_angle_alpha" ${folder}/test.cif | awk '{print $2}')
  beta=$(grep "_cell_angle_beta" ${folder}/test.cif | awk '{print $2}')
  gamma=$(grep "_cell_angle_gamma" ${folder}/test.cif | awk '{print $2}')
  cell_parameters=$(echo "$a $b $c $alpha $beta $gamma")
  echo $scale $energy $gnorm ${cell_parameters} >> cell_energy_expansion_gulp.txt
 fi
done
sort -gk1 cell_energy_expansion_gulp.txt | sed 's/CRYST1 //g' > c
mv c cell_energy_expansion_gulp.txt
exit 0
