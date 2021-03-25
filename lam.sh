#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

rm res/resultat.txt

liste_fichiers=`ls 20x20`
for fichier in $liste_fichiers
do
chaine="20x20/"
echo $chaine$fichier
mpirun -np 4 ./bin/AATESTPROG --namelist fort.4.20x20 --field-file $chaine$fichier --time 1 > AATESTPROG.eo 2>1
cp AATESTPROG.eo  res/AATESTPROG.eo$fichier
cd res
grep '0.000000000000000' AATESTPROG.eo$fichier | wc -l >> resultat.txt
cd ..
done

#mpirun -np 1 ./bin/AATESTPROG --namelist fort.4.20x20 --field-file 20x20/AATESTPROG.20x20.gp.000045.dat --time 1 > AATESTPROG.eo 2>&1
