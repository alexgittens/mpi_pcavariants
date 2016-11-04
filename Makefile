#SHELL="/bin/bash"
.SHELLARGS="-l -c"

# All parameters chosen for Cori

# use 7 slurm nodes
MEDIUMFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/ocean.h5
MEDIUMVARNAME=rows
MEDIUMLATFNAME=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/observedLatitudes.csv
LARGEMETADATADIR=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/
MEDIUMNUMROWS=6177583
MEDIUMNUMCOLS=8096
MEDIUMRANK=100
MEDIUMOUTFNAME=mediumeofs.h5
MEDIUMEOF3DFNAME=mediumoefs.nc

# use 40 slurm nodes : takes less than 12 minutes
# NB: the determining factor in the number of nodes is that Parallel-HDF5 can only load less than 2GB in each call to h5fread
LARGEFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/output/ocean.h5
LARGEVARNAME=rows
LARGELATFNAME=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/output/observedLatitudes.csv
LARGEMETADATADIR=/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/output/
LARGENUMROWS=6177583
LARGENUMCOLS=47112
LARGERANK=100
LARGEOUTFNAME=largeeofs.h5
LARGEEOF3DFNAME=largeeofs.nc

# use 30 slurm nodes
# NB: this has to hold both temp and rho in memory
CESMFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/CESM_conversion/output/cesm.h5
CESMTEMPDATASETNAME=temp
CESMRHODATASETNAME=rho
CESMLATFNAME=/global/cscratch1/sd/gittens/conversion-code/CESM_conversion/output/observedLatitudes.csv
CESMMETADATADIR=/global/cscratch1/sd/gittens/conversion-code/CESM_conversion/output/
CESMRANK=100
CESMOUTFNAME=cesmeofs.h5
CESMEOF3DFNAME=cesmeofs.nc

# use 30 CORI nodes
THERMOCLINEFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/thermoclineOcean.h5
THERMOCLINEDATASETNAME=rows
THERMOCLINELATFNAME=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/observedLatitudes.csv
THERMOCLINEMETADATADIR=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/
THERMOCLINERANK=100
THERMOCLINEOUTFNAME=thermoclineeofs.h5
THERMOCLINEEOF3DFNAME=thermoclineeofs.nc
all: cori

cori: pca.c
	module load cray-hdf5-parallel && \
	cc -g -o pca pca.c -larpack -L.

coritest: cori
	module load cray-hdf5-parallel && \
	srun -u -n 5 ./pca test2.hdf5 temperatures rowWeights.csv 10 4 3 out.hdf5
	
mediumscale: cori
	module load cray-hdf5-parallel && \
	srun -u -c 1 -n 187 ./pca ${MEDIUMFNAMEIN} ${MEDIUMVARNAME} ${MEDIUMLATFNAME} ${MEDIUMNUMROWS} ${MEDIUMNUMCOLS} ${MEDIUMRANK} ${MEDIUMOUTFNAME} && \
	python reshape.py ${MEDIUMOUTFNAME} "latweighting+centering" ${MEDIUMMETADATADIR} ${MEDIUMEOF3DFNAME}

cfsro: cori
	module load python cray-hdf5-parallel && \
	srun -u -c 1 -n 1252 ./pca ${LARGEFNAMEIN} ${LARGEVARNAME} ${LARGELATFNAME} ${LARGENUMROWS} ${LARGENUMCOLS} ${LARGERANK} ${LARGEOUTFNAME} && \
	python cfsroreshape.py ${LARGEOUTFNAME} "latweighting+centering" ${LARGEMETADATADIR} ${LARGEEOF3DFNAME}

# deprecated, but keep around in case need to compile on Edison at some point
edison: pca.c
	module load cray-hdf5-parallel && \
	cc -o pca pca.c -I$$HDF5_INCLUDE_OPTS -static -lhdf5 -larpack -L. -L$$CRAY_LD_LIBRARY_PATH

modular: computations.c io.c pca.h modularpca.c 
	module load cray-hdf5-parallel && \
	cc -std=c99 -g -o modularpca modularpca.c computations.c io.c -larpack -I. -L.

cesm: modular
	module load python cray-hdf5-parallel && \
	srun -u -c 1 -n 960 ./modularpca ${CESMFNAMEIN} ${CESMFNAMEIN} ${CESMTEMPDATASETNAME} ${CESMRHODATASETNAME} ${CESMLATFNAME} ${CESMRANK} ${CESMOUTFNAME} && \
	python cesmreshape.py ${CESMOUTFNAME} "latweighting+centering" ${CESMMETADATADIR} ${CESMEOF3DFNAME}

thermocline: modular
	module load python cray-hdf5-parallel && \
	srun -u -c 1 -n 960 ./modularpca  ${THERMOCLINEFNAMEIN} NULL ${THERMOCLINEDATASETNAME} NULL ${THERMOCLINELATFNAME} ${THERMOCLINERANK} ${THERMOCLINEOUTFNAME} && \
	python thermoclinereshape.py ${THERMOCLINEOUTFNAME} "latweighting+centering" ${THERMOCLINEMETADATADIR} ${THERMOCLINEEOF3DFNAME}
