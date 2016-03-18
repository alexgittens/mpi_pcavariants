#include "hdf5.h"
#include "mpi.h"
#include <stdio.h>
#include "cblas.h"

//#define INFILENAME "/global/cscratch1/sd/gittens/CFSROhdf5/oceanTemps.hdf5"
#define INFILENAME "test.hdf5"
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
    //int numcols = 46715;
    //int numrows = 6349676;
    int numcols = 4;
    int numrows = 10;

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

    // compute C = A'*B = A' * (A*Omega)
    // compute C = A'*B
    // a la https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/280875
    /*
    int height = 2;
    int numcols = 4;
    int r = 3;
    double * Alocal, * Scratch, * Clocal, * Omega;
    Omega = (double *) malloc(numcols*r*sizeof(double));
    Alocal = (double *) malloc(height*numcols*sizeof(double));
    Scratch = (double *) malloc(height*r*sizeof(double));
    Clocal = (double *) malloc(numcols*r*sizeof(double));

    int rowIdx;
    int colIdx;
    for (rowIdx = 0; rowIdx < numcols; rowIdx = rowIdx + 1) { 
        for(colIdx = 0; colIdx < r; colIdx = colIdx + 1) {
            Omega[rowIdx*r + colIdx] = rowIdx + colIdx/2.0;
            printf("Omega[%d][%d] = %f\n", rowIdx, colIdx, Omega[rowIdx*r + colIdx]);
        }
    }

    for (rowIdx = 0; rowIdx < height; rowIdx = rowIdx + 1) { 
        for(colIdx = 0; colIdx < numcols; colIdx = colIdx + 1) {
            if (colIdx == 2 && rowIdx == 1) {
                Alocal[rowIdx*numcols + colIdx] = 5; 
            } else {
                Alocal[rowIdx*numcols + colIdx] = rowIdx + 1;
            }
            printf("A[%d][%d] = %f\n", rowIdx, colIdx, Alocal[rowIdx*numcols + colIdx]);
        }
    }
    multiplyGramianChunk(Alocal, Omega, Clocal, Scratch, height, numcols, r);

    if(mpi_rank == 0) {
        for (rowIdx = 0; rowIdx < numcols; rowIdx = rowIdx + 1) { 
            for(colIdx = 0; colIdx < r; colIdx = colIdx + 1) {
                printf("C[%d][%d] = %f\n", rowIdx, colIdx, Clocal[rowIdx*r + colIdx]);
            }
        }
    }

    */

    /* Open the file */
    file_id = H5Fopen(INFILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Load my portion of the data */
    double * Alocal = (double *) malloc( mynumrows * numcols * sizeof(double));

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

    /*
    int rowIdx, colIdx;

    for (rowIdx = 0; rowIdx < mynumrows; rowIdx = rowIdx + 1) {
        for (colIdx = 0; colIdx < numcols; colIdx = colIdx + 1) {
            printf("A[%d][%d] = %f \n", rowIdx + startingrow, colIdx, Alocal[rowIdx*numcols + colIdx]);
        }
    }
    */

    H5Dclose(dataset_id);
    H5Fclose(file_id);

    MPI_Finalize();
    return 0;
}

/* computes A*(A'*Omega) = A*S , so Scratch must have size rowsA*colsOmega */
void multiplyGramianChunk(double A[], double Omega[], double C[], double Scratch[], int rowsA, int colsA, int colsOmega) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowsA, colsOmega, colsA, 1.0, A, colsA, Omega, colsOmega, 0.0, Scratch, colsOmega);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, colsA, colsOmega, rowsA, 1.0, A, colsA, Scratch, colsOmega, 0.0, C, colsOmega);
}
