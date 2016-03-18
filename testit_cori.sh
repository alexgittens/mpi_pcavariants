#!/bin/bash -l
module load hdf5
cc testcase.c $HDF5_INCLUDE_OPTS -static -lhdf5 $HDF5_POST_LINK_OPTS
srun ./a.out
