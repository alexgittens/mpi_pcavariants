#include "hdf5.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cblas.h"
#include "lapacke.h"
#include <math.h>
#include "pca.h"
#include <limits.h>
#include <stdbool.h>

// NB: CBLAS has nonconstant overhead, because after operations, it stores the output in row major
// TODO : use BLAS level 2 to compute the matrix-vector products!
// TODO : allow each process to load more than 2GB at a time by doing reading in chunks (will make more scalable: really, the whole point of doing 
// distributed matrix-vector products is to allow the matrix to be fit in memory, so we should make ourselves memory-limited)
// TODO : use memory-aligned mallocs

/*
#define MAX_MATVECPRODS 10000000
#define MAX_RESTARTS 100000
*/
#define MAX_MATVECPRODS INT_MAX
#define MAX_RESTARTS INT_MAX
#define CONVERGENCE_TOL 1e-13

int main(int argc, char **argv) {

    /* MPI variables */
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    /* Initialize MPI */
    int mpi_size, mpi_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    double totalTime = MPI_Wtime();

    char * temp_infilename = argv[1];
    char * rho_infilename = argv[2];
    char * temp_datasetname = argv[3];
    char * rho_datasetname = argv[4];
    char * weightsFname = argv[5];
    int numeigs = atoi(argv[6]);
    char * outfname = argv[7];

    if (mpi_rank == 0) {
        printf("Arguments:\n");
        printf("\tInput temperatures filename: %s\n", temp_infilename);
        printf("\tInput density filename: %s\n", rho_infilename);
        printf("\tInput temperatures dataset: %s\n", temp_datasetname);
        printf("\tInput density dataset: %s\n", rho_datasetname);
        printf("\tRow to latitude mapping file: %s\n", weightsFname);
        printf("\tTarget rank: %d\n", numeigs);
        printf("\tOutput filename: %s\n", outfname);
    }

    /* Load the temperature and density datasets from file */
    distMatrixInfo *matInfo = (distMatrixInfo *) malloc( sizeof(distMatrixInfo));
    distGatherInfo *eigInfo = (distGatherInfo *) malloc( sizeof(distGatherInfo));
    double *localRowChunk, *localRhoRowChunk;
    double readDimsTime = getMatrixInfo(temp_infilename, temp_datasetname, &comm, &info, matInfo);
    double readTempTime = loadMatrix(temp_infilename, temp_datasetname, matInfo, &localRowChunk, &comm, &info);
    double readRhoTime = 0;
    if (strcmp(rho_datasetname, "NULL")) {
        readRhoTime = loadMatrix(rho_infilename, rho_datasetname, matInfo, &localRhoRowChunk, &comm, &info);
    }

    /* Load the row weights (sqrt of cos of the latitudes corresponding to each row) from file */
    double *rowWeights = (double *) malloc( sizeof(double)*matInfo->numrows );
    double readWeightsTime;
    if ( mpi_rank == 0) {
        readWeightsTime = loadRowWeights(weightsFname, rowWeights);
    }
    MPI_Bcast(rowWeights, matInfo->numrows, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&readWeightsTime, 1, MPI_DOUBLE, 0, comm);
    if (mpi_rank == 0) {
        printf("Time to read input datasets: %f\n", readDimsTime + readTempTime + readRhoTime + readWeightsTime);
    }

    genGatherInfo(matInfo, numeigs, eigInfo);
    int numcols = matInfo->numcols;
    int numrows = matInfo->numrows;

    /* Compute row means and center the rows of the temperature data:
     * first mutiply the temperature by density, then compute and subtract mean temperatures, then rescale by latitude weights
     */
    double *meanVec = malloc(numrows * sizeof(double));
    double preprocessingTime = MPI_Wtime();
    if (strcmp(rho_datasetname, "NULL")) {
        dhad(localRowChunk, localRhoRowChunk, matInfo->localrows * numcols); //hadamard product: for CESM, multiply temperature by density
    }
    computeAndSubtractRowMeans(localRowChunk, meanVec, matInfo);
    rescaleRows(localRowChunk, rowWeights, matInfo);

    preprocessingTime = MPI_Wtime() - preprocessingTime;
    if (mpi_rank == 0) {
        printf("Time to preprocess the data (center rows and weight rows by sqrt(cos(lat))): %f\n", preprocessingTime);
    }

    double frobnorm, tcompfrobstr, tcompfrobstp;
    tcompfrobstr = MPI_Wtime();
    computeFrobNorm(localRowChunk, &frobnorm, matInfo);
    tcompfrobstp = MPI_Wtime();
    if (mpi_rank == 0) {
        printf("Time to compute Frobnenius norm of A: %f\n", tcompfrobstp - tcompfrobstr);
        printf("Frobenius norm of A: %g\n", frobnorm);
    }

    /* Allocate space for the vector used in Lanczos iterations and the (local to this process) scratch matrices used 
     * when computing the Gram-vector product */
    double * vector = (double *) malloc( numcols * sizeof(double));
    scratchMatrices * scratchSpace = (scratchMatrices *) malloc( sizeof(scratchMatrices));
    scratchSpace->Scratch = (double *) malloc( matInfo->localrows * sizeof(double));
    scratchSpace->Scratch2 = (double *) malloc( numcols * sizeof(double));
    scratchSpace->Scratch3 = (double *) malloc( matInfo->localrows * numeigs * sizeof(double));
    double * singVals = (double *) malloc( numeigs * sizeof(double));
    double * rightSingVecs = (double *) malloc( numeigs * numcols * sizeof(double));
    if (vector == NULL || scratchSpace == NULL || scratchSpace->Scratch == NULL || 
        scratchSpace->Scratch2 == NULL || scratchSpace->Scratch3 == NULL || 
        singVals == NULL || rightSingVecs == NULL) {
        printf("Out of memory on process %d\n", mpi_rank);
        exit(-1);
    }

    if(mpi_rank == 0) {
        printf("Computing the EVD of the Gram matrix\n");
    }

    /* Define ARPACK working variables and parameters*/
    int ido = 0;
    int ncv = 2*numeigs > numcols ? numcols : 2*numeigs; // ncv > nev and ncv < n (but ncv >= 2*nev recommended)
    double tol = CONVERGENCE_TOL;
    double * resid = (double *) malloc( numcols * sizeof(double));
    double * v = (double *) malloc(numcols * ncv *sizeof(double));
    int iparam[11] = {1, 0, 30, 1, 0, 0, 1, 0, 0, 0, 0};
    iparam[2] = MAX_RESTARTS;
    int ipntr[11];
    double * workd = (double *) malloc(3*numcols*sizeof(double));
    int lworkl = ncv*(ncv + 8);
    double * workl = (double *) malloc(lworkl*sizeof(double));
    int arpack_info = 0;
    char bmat = 'I';
    char which[3] = "LM";
    if ( resid == NULL || v == NULL || workd == NULL || workl == NULL) {
        printf("Out of memory on process %d\n", mpi_rank);
        exit(-1);
    }

    /* initialize ARPACK */
    if (mpi_rank == 0) {
        dsaupd_(&ido, &bmat, &numcols, which,
                &numeigs, &tol, resid,
                &ncv, v, &numcols,
                iparam, ipntr, workd,
                workl, &lworkl, &arpack_info);
        cblas_dcopy(numcols, workd + ipntr[0] - 1, 1, vector, 1); 
    }
    MPI_Bcast(&ido, 1, MPI_INTEGER, 0, comm);
    MPI_Bcast(vector, numcols, MPI_DOUBLE, 0, comm);

    /* Call ARPACK until convergence */
    int numMatVecProds = 0;
	double tgrammv = 0., grammvstr, grammvstp;
	double tarpk = 0., arpkstr, arpkstp;
    bool continueFlag = true;
	while(numMatVecProds < MAX_MATVECPRODS && continueFlag) {
        if(numMatVecProds % 100 == 0 && mpi_rank == 0) {
            printf("Computed %d matrix-vector products against the grammian\n", numMatVecProds);
        }
        if (ido == 1 || ido == -1) {
            grammvstr = MPI_Wtime();
            distributedGramianVecProd(localRowChunk, vector, matInfo, scratchSpace);
            grammvstp = MPI_Wtime();
            tgrammv += grammvstp - grammvstr;
            arpkstr = MPI_Wtime();
            if (mpi_rank == 0) {
                cblas_dcopy(numcols, vector, 1, workd + ipntr[1] - 1, 1); // y = A x
                dsaupd_(&ido, &bmat, &numcols, which,
                        &numeigs, &tol, resid,
                        &ncv, v, &numcols,
                        iparam, ipntr, workd,
                        workl, &lworkl, &arpack_info);
                cblas_dcopy(numcols, workd + ipntr[0] - 1, 1, vector, 1);
                if (iparam[4] < numeigs) {
                    continueFlag = true;
                } else {
                    continueFlag = false;
                }
            }
            MPI_Bcast(vector, numcols, MPI_DOUBLE, 0, comm); 
            MPI_Bcast(&continueFlag, 1, MPI_C_BOOL, 0, comm);
        }
        MPI_Bcast(&ido, 1, MPI_INTEGER, 0, comm);
        arpkstp = MPI_Wtime();
        tarpk += arpkstp - arpkstr;
        numMatVecProds++;
	}

    if (mpi_rank == 0) {
        if (numMatVecProds == MAX_MATVECPRODS) {
            printf("Terminated eval decomposition due to reaching max_iters mat-vec products: %d\n", MAX_MATVECPRODS);
        } else {
            printf("Completed eval decomposition after %d mat-vec products\n", numMatVecProds);
        }
    }

	if(mpi_rank == 0){
		printf("Time to perform distributed Gram matrix-vectors: %f\n", tgrammv); 
		printf("Time spend in arpack: %f\n", tarpk); 
	}


    if (mpi_rank == 0) {
        int numTotalMatVecProds = iparam[8];
        int numEigsComputed = iparam[4];
        double trtzstr = MPI_Wtime();
		printf("Used %d matrix-vector products to converge to %d eigenvalue\n", numTotalMatVecProds, numEigsComputed);
        printf("Extracting the right singular vectors\n");

        int rvec = 1; // compute Ritz vectors
        char HowMny = 'A';
        int * select = (int * ) malloc(ncv * sizeof(int));
        double sigma = 0;
    	
        // eigenvalues and eigenvectors are returned in ascending order
        // eigenvectors are returned in column major form
        double * svtranspose = (double *) malloc( numeigs * numcols * sizeof(double));
		if(svtranspose == NULL || select == NULL)
        	printf("Uh Oh, svtranspose is NULL\n");
			
		dseupd_(&rvec, &HowMny, select,
                singVals, svtranspose, &numcols,
                &sigma, &bmat, &numcols,
                which, &numeigs, &tol, 
                resid, &ncv, v,
                &numcols, iparam, ipntr,
                workd, workl, 
                &lworkl, &arpack_info);

        mattrans(svtranspose, numeigs, numcols, rightSingVecs);
        flipcolslr(rightSingVecs, numcols, numeigs);

		double trtzstp = MPI_Wtime();
		printf("Time to compute Ritz vectors: %f\n", trtzstp - trtzstr);
        free(svtranspose); 
        free(select);
    }

    if (mpi_rank ==0) {
        printf("Computing AV\n");
    }
    double tcompavstr, tcompavstp;
	tcompavstr = MPI_Wtime();
	MPI_Bcast(rightSingVecs, numeigs*numcols, MPI_DOUBLE, 0, comm);
    double * AV = (double *) malloc( numrows * numeigs * sizeof(double));
    distributedMatMatProd(localRowChunk, rightSingVecs, AV, matInfo, eigInfo, scratchSpace);
	tcompavstp = MPI_Wtime();
    if (mpi_rank == 0) {
		printf("Time to compute AV: %f\n", tcompavstp - tcompavstr);
    }

    // dgesdd returns its singular values in descending order
    // note that the right singular vectors of AV should by definition be the identity 
    if (mpi_rank == 0) {

        printf("Computing the SVD of AV\n");
        double svdTime = MPI_Wtime();
		double * U = (double *) malloc( numrows * numeigs * sizeof(double));
        double * VT = (double *) malloc( numeigs * numeigs * sizeof(double));
        double * V = (double *) malloc( numeigs * numeigs * sizeof(double));
        double * finalV = (double *) malloc( numcols * numeigs * sizeof(double));
        double * singvals = (double *) malloc( numeigs * sizeof(double));
        LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', numrows, numeigs, AV, numeigs, singvals, U, numeigs, VT, numeigs);
        mattrans(VT, numeigs, numeigs, V);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numcols, numeigs, numeigs, 1.0, rightSingVecs, numeigs, V, numeigs, 0.0, finalV, numeigs);
        svdTime = MPI_Wtime() - svdTime;
		printf("Time to compute SVD of AV: %f\n", svdTime);

        // Write the output
        printf("Writing the output of the SVD to file\n");
        writeSVD(outfname, eigInfo, U, finalV, singvals, meanVec, rowWeights);

        free(U);
        free(VT);
        free(V);
        free(finalV);
        free(singvals);
    }

    free(localRowChunk);
    if (strcmp(rho_datasetname, "NULL")) {
        free(localRhoRowChunk);
    }
    free(rowWeights);
    free(vector);
    free(resid);
    free(v);
    free(workd);
    free(workl);
    freeMatrixInfo(matInfo);
    freeScratchMatrices(scratchSpace);
    if (mpi_rank == 0) {
        freeGatherInfo(eigInfo);
    }

    totalTime = MPI_Wtime() - totalTime;
	if(mpi_rank == 0) {
		printf("Total PCA elapsed time (including Frob norm calculation): %f\n", totalTime);
    }
	MPI_Finalize();
    return 0;
}
