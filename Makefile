#SHELL="/bin/bash"
.SHELLARGS="-l -c"

# All parameters chosen for Cori (128GB/node, 32 cores/node), and keeping in mind that parallel-hdf5 can only load less than 2GB in each call to h5fread

# use 40 slurm nodes
CFSROFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/output/ocean.h5
CFSRODATASETNAME=rows
CFSROLATFNAME=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/output/observedLatitudes.csv
CFSROMETADATAFNAME=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/output/oceanMetadata.npz
CFSRORANK=100
CFSROOUTFNAME=cfsroeofs.h5
CFSROEOF3DFNAME=cfsroeofs.nc

# use 30 CORI nodes
THERMOCLINEFNAMEIN=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/thermoclineOcean.h5
THERMOCLINEDATASETNAME=rows
THERMOCLINELATFNAME=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/observedLatitudes.csv
THERMOCLINEMETADATAFNAME=/global/cscratch1/sd/gittens/conversion-code/CFSRO_conversion/thermoclineOutput/thermoclineOceanMetadata.npz
THERMOCLINERANK=100
THERMOCLINEOUTFNAME=thermoclineeofs.h5
THERMOCLINEEOF3DFNAME=thermoclineeofs.nc

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

all: modular

modular: computations.c io.c pca.h modularpca.c 
	module load cray-hdf5-parallel && \
	cc -std=c99 -g -o modularpca modularpca.c computations.c io.c -larpack -lm -I. -L. 

cfsro: modular
	#module load python cray-hdf5-parallel && \
	#srun -u -c 1 -n 1252 ./modularpca ${CFSROFNAMEIN} NULL ${CFSRODATASETNAME} NULL ${CFSROLATFNAME} ${CFSRORANK} ${CFSROOUTFNAME} && \
	python cfsroreshape.py ${CFSROOUTFNAME} "latweighting+centering" ${CFSROMETADATAFNAME} ${CFSROEOF3DFNAME}

thermocline: modular
	module load python cray-hdf5-parallel && \
	srun -u -c 1 -n 960 ./modularpca  ${THERMOCLINEFNAMEIN} NULL ${THERMOCLINEDATASETNAME} NULL ${THERMOCLINELATFNAME} ${THERMOCLINERANK} ${THERMOCLINEOUTFNAME} && \
	python thermoclinereshape.py ${THERMOCLINEOUTFNAME} "latweighting+centering" ${THERMOCLINEMETADATAFNAME} ${THERMOCLINEEOF3DFNAME}

cesm: modular
	module load python cray-hdf5-parallel && \
	srun -u -c 1 -n 960 ./modularpca ${CESMFNAMEIN} ${CESMFNAMEIN} ${CESMTEMPDATASETNAME} ${CESMRHODATASETNAME} ${CESMLATFNAME} ${CESMRANK} ${CESMOUTFNAME} && \
	python cesmreshape.py ${CESMOUTFNAME} "latweighting+centering" ${CESMMETADATADIR} ${CESMEOF3DFNAME}

# deprecated, but keep around in case need to compile on Edison at some point
edison: pca.c
	module load cray-hdf5-parallel && \
	cc -o pca pca.c -I$$HDF5_INCLUDE_OPTS -static -lhdf5 -larpack -L. -L$$CRAY_LD_LIBRARY_PATH
