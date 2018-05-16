#!/bin/bash
if [ -f cell_energy_expansion_gulp.txt ] ; then rm cell_energy_expansion_gulp.txt ; touch cell_energy_expansion_gulp.txt ; fi
temperature=50.0
pressure=0.0
nCPU=4
while read line ; do
 file=$(echo $line | awk '{print $1}')
 energy_VASP=$(echo $line | awk '{print $2}')
 structure=$(echo $file | sed 's/struc\///g' | sed 's/\.cif//g' )
 folder=${structure}_dir
 if [ ! -d $folder ] ; then
  mkdir $folder
  cd $folder
   cp ../*.f90 ../peros.lib .
   gfortran cif2lammps.f90 -o cif2lammps
   cp ../../struc/${structure}.cif .
   ./cif2lammps -R -S -c ${structure}.cif > ${structure}.l2gout
   mpirun --np $nCPU ~/bin/gulp < ${structure}.gin > ${structure}.gout 
  cd ..
 else 
  energy=$(grep "Total lattice energy       =" ${folder}/${structure}.gout | grep "eV" | awk '{print $5}' | tail -n1 )
  gnorm=0.0
  a=$(grep "_cell_length_a" ${folder}/test.cif | awk '{print $2}')
  b=$(grep "_cell_length_b" ${folder}/test.cif | awk '{print $2}')
  c=$(grep "_cell_length_c" ${folder}/test.cif | awk '{print $2}')
  alpha=$(grep "_cell_angle_alpha" ${folder}/test.cif | awk '{print $2}')
  beta=$(grep "_cell_angle_beta" ${folder}/test.cif | awk '{print $2}')
  gamma=$(grep "_cell_angle_gamma" ${folder}/test.cif | awk '{print $2}')
  cell_parameters=$(echo "$a $b $c $alpha $beta $gamma")
  echo ${energy_VASP} $energy $gnorm ${cell_parameters} >> cell_energy_expansion_gulp.txt
 fi
done < list
sort -gk1 cell_energy_expansion_gulp.txt | sed 's/CRYST1 //g' > c
mv c cell_energy_expansion_gulp.txt
exit 0
