#!/bin/bash -l
module load cray-hdf5
cc testcase.c -I$HDF5_INCLUDE_OPTS -static -lhdf5 -L$CRAY_LD_LIBRARY_PATH
srun -n 5 ./a.out
