SHELL="/bin/bash"
.SHELLARGS="-l -c"
FNAMEIN="/global/cscratch1/sd/jialin/climate/oceanTemps.hdf5"

all: cori

cori: pca.c
	module load hdf5-parallel; \
	echo "$$HDF5_INCLUDE_OPTS"; \
	echo "$$HDF5_POST_LINK_OPTS"; \
	cc -o pca pca.c "$$HDF5_INCLUDE_OPTS" -static -lhdf5 -larpack -L. "$$HDF5_POST_LINK_OPTS"

edison: pca.c
	module load cray-hdf5-parallel; \
	cc -o pca pca.c -I$$HDF5_INCLUDE_OPTS -static -lhdf5 -larpack -L. -L$$CRAY_LD_LIBRARY_PATH

coritest: cori
	modul load hdf5-parallel; \
	srun -u -n 5 ./pca test2.hdf5 temperatures 10 4 3 out.hdf5

edisontest: edison
	srun -u -n 5 ./pca test2.hdf5 temperatures 10 4 3 out.hdf5
	
coriproduction: cori
	module load hdf5-parallel; \
	srun -u -c 1 -n 60 --mem-per-cpu=51200 ./pca ${FNAMEIN} temperatures 6349676 46715 20 out.hdf5
