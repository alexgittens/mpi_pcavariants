#SHELL="/bin/bash"
.SHELLARGS="-l -c"

# use 7 slurm nodes
MEDIUMFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/ocean.h5
MEDIUMVARNAME=rows
MEDIUMLATFNAME=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/observedLatitudes.csv
MEDIUMNUMROWS=6177583
MEDIUMNUMCOLS=8096
MEDIUMRANK=100
MEDIUMOUTFNAME=mediumeofs.h5

# use 40 slurm nodes
LARGEFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/output/ocean.h5
LARGEVARNAME=rows
LARGELATFNAME=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/output/observedLatitudes.csv
LARGENUMROWS=6177583
LARGENUMCOLS=47112
LARGERANK=100
LARGEOUTFNAME=largeeofs.h5

all: cori

cori: pca.c
	module load cray-hdf5-parallel; \
	cc -g -o pca pca.c -larpack -L.

edison: pca.c
	module load cray-hdf5-parallel; \
	cc -o pca pca.c -I$$HDF5_INCLUDE_OPTS -static -lhdf5 -larpack -L. -L$$CRAY_LD_LIBRARY_PATH

coritest: cori
	module load cray-hdf5-parallel; \
	srun -u -n 5 ./pca test2.hdf5 temperatures rowWeights.csv 10 4 3 out.hdf5

edisontest: edison
	srun -u -n 5 ./pca test2.hdf5 temperatures 10 4 3 out.hdf5
	
mediumscale: cori
	module load cray-hdf5-parallel; \
	srun -u -c 1 -n 187 ./pca ${MEDIUMFNAMEIN} ${MEDIUMVARNAME} ${MEDIUMLATFNAME} ${MEDIUMNUMROWS} ${MEDIUMNUMCOLS} ${MEDIUMRANK} ${MEDIUMOUTFNAME}

largescale: cori
	module load cray-hdf5-parallel; \
	srun -u -c 1 -n 1252 ./pca ${LARGEFNAMEIN} ${LARGEVARNAME} ${LARGELATFNAME} ${LARGENUMROWS} ${LARGENUMCOLS} ${LARGERANK} ${LARGEOUTFNAME}
