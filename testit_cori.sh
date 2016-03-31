#!/bin/bash -l
module load hdf5-parallel
cc pca.c $HDF5_INCLUDE_OPTS -static -lhdf5 -larpack -L. $HDF5_POST_LINK_OPTS
srun -u -n 5 ./a.out test2.hdf5 temperatures 10 4 3 out.hdf5
