#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

mpirun -np 2 ./bin/AATESTPROG --namelist fort.4.t31 --time 1 > AATESTPROG.eo 2>&1

/home/ms/fr/sor/3d/glgrib/glgrib.sh AATESTPROG.fa%ZZZFFFFF
mv snapshot_0000.png ~/tmp/.

ssh ecgate scp ~/tmp/snapshot_0000.png phi001@90.76.140.145:tmp/snapshot_0000.png



