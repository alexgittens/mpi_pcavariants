#include "hdf5.h"
#include "mpi.h"
#include <stdio.h>
#include "cblas.h"

#define INFILENAME "/global/cscratch1/sd/gittens/CFSROhdf5/oceanTemps.hdf5"
//#define INFILENAME "test.hdf5"
#define DATASET "temperatures"

void multiplyGramianChunk(double A[], double Omega[], double C[], double Scratch[], int rowsA, int colsA, int colsOmega); 

// optimizations: use memory-aligned mallocs

int main(int argc, char **argv) {

    /* HDF5 API definitions */
    hid_t file_id, dataset_id, dataspace, memspace; 
    herr_t status; 
    hsize_t offset[2], count[2], offset_out[2];

    /* MPI variables */
    int mpi_size, mpi_rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    /* Allocate the correct portion of the input to each processor */
    int numcols = 46715;
    int numrows = 6349676;
    //int numcols = 4;
    //int numrows = 10;

    int rowsModSize = numrows/mpi_size;
    int numWithMoreRows = numrows % mpi_size;
    int startingrow, mynumrows;

    if (mpi_rank < numWithMoreRows) {
        startingrow = mpi_rank*(rowsModSize + 1);
        mynumrows = rowsModSize + 1;
    } else {
        startingrow = numWithMoreRows*(rowsModSize + 1) + (mpi_rank - numWithMoreRows)*rowsModSize;
        mynumrows = rowsModSize;
    }
    printf("Rank %d: assigned %d rows, %d--%d\n", mpi_rank, mynumrows, startingrow, startingrow + mynumrows - 1);

    /* Load my portion of the data */
    double * Alocal = (double *) malloc( mynumrows * numcols * sizeof(double));

    // Takes about 12 minutes to load the dataset
    int rem;
    int ioparallelism = 50;
    for (rem = 0; rem < ioparallelism; rem = rem + 1) {
        if (mpi_rank % ioparallelism == rem) {
            // Open the file 
            file_id = H5Fopen(INFILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

            dataset_id = H5Dopen2(file_id, DATASET, H5P_DEFAULT);
            dataspace = H5Dget_space(dataset_id);
            offset[0] = startingrow;
            offset[1] = 0;
            count[0] = mynumrows;
            count[1] = numcols;
            status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

            offset_out[0] = 0;
            offset_out[1] = 0;
            memspace = H5Screate_simple(2, count, NULL);
            status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL);
            status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, Alocal);

            printf("Rank %d: loaded my data\n", mpi_rank);

            H5Dclose(dataset_id);
            H5Fclose(file_id);
        }
        MPI_Barrier(comm);
     }

    MPI_Finalize();
    return 0;
}

/* computes A*(A'*Omega) = A*S , so Scratch must have size rowsA*colsOmega */
void multiplyGramianChunk(double A[], double Omega[], double C[], double Scratch[], int rowsA, int colsA, int colsOmega) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowsA, colsOmega, colsA, 1.0, A, colsA, Omega, colsOmega, 0.0, Scratch, colsOmega);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colsA, colsOmega, rowsA, 1.0, A, colsA, Scratch, colsOmega, 0.0, C, colsOmega);
}
