/*
**  PROGRAM: Matrix Multiply
**
**  PURPOSE: This is a simple matrix multiply program. 
**           It will compute the product
**
**                C  = A * B
**
**           A and B are set to constant matrices so we
**           can make a quick test of the multiplication.
**
**  HISTORY: Written by Tim Mattson, Nov 1999.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define AVAL 3.0
#define BVAL 5.0
#define TOL  0.001

int main(int argc, char **argv)
{
	int Ndim, Pdim, Mdim;   /* A[N][P], B[P][M], C[N][M] */
	int i,j,k;
	double *A, *B, *C, cval, err, errsq;
      double dN, mflops;
	double start_time, run_time;
    int rank, size;

    if (argc < 2) {
        printf("Usage: %s <matrix_order>\n", argv[0]);
        return -1;
    }

	int order = atoi(argv[1]);
	Ndim = order;
	Pdim = order;
	Mdim = order;

	/* Initialize matrices */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    B = (double *)malloc(Pdim*Mdim*sizeof(double));
    if (rank == 0) {
        A = (double *)malloc(Ndim*Pdim*sizeof(double));
        C = (double *)malloc(Ndim*Mdim*sizeof(double));

        for (i=0; i<Ndim; i++)
            for (j=0; j<Pdim; j++)
                *(A+(i*Ndim+j)) = AVAL;

        for (i=0; i<Pdim; i++)
            for (j=0; j<Mdim; j++)
                *(B+(i*Pdim+j)) = BVAL;

        for (i=0; i<Ndim; i++)
            for (j=0; j<Mdim; j++)
                *(C+(i*Ndim+j)) = 0.0;
    }

    start_time = MPI_Wtime();
    /* =========================== AREA OF INTEREST ============================ */

    int root = 0;
    MPI_Comm world_comm = MPI_COMM_WORLD;
    int A_rows = Ndim;
    int A_cols = Pdim;

    int B_rows = Pdim;
    int B_cols = Mdim;

    int C_rows = Ndim;
    int C_cols = Mdim;

    // send matrix B to all ranks
    int B_values = B_rows * B_cols;
    MPI_Bcast(B, B_values, MPI_DOUBLE, root, world_comm);

    // scatter A rows to all ranks
    int total_ranks = size;
    int local_rows =  A_rows / total_ranks; 
    int local_A_values = local_rows * A_cols; //numero di elementi di A, quindi numero di righe * Acols

    double* local_A = (double*) malloc(local_A_values * sizeof(double));
    MPI_Scatter(A, local_A_values, MPI_DOUBLE, local_A, local_A_values, MPI_DOUBLE, root, world_comm);


    //alloca C = local_rows * C_cols e compute locale
    //approccio letterale di moltiplicazione
    int local_C_values = local_rows * C_cols;
    double* local_C = (double*) calloc(local_C_values, sizeof(double)); //clean mem

    // for (int row = 0; row < local_rows; row++){ // righe di a locale
    //     for (int col = 0; col < B_cols; col ++){ // colonne di b, c
    //         for (int k = 0; k < A_cols; k++){   //somma lungo la riga di A + colonna di B
    //             local_C[row * C_cols + col] += local_A[row * A_cols + k] * B[k * B_cols + col];
    //             //ps. è sempre colonne totali * riga attaule + col attuale
    //         }
    //     }
    // }


    //ottimizzato
    for (int i = 0; i < local_rows; i++){
        for(int j = 0; j < A_cols; j++){ //itera su righe di B
            double a = local_A[A_cols * i + j];
            for (int k = 0; k < B_cols; k ++){ //itera su col di B / C
                double b = B[B_cols * j + k];
                local_C[C_cols * i + k] += a * b; //c ha la stessa riga di a
            }
        }
    }




    /* Do the matrix product */
    // for (i=0; i<Ndim; i++){
	// 	for (j=0; j<Mdim; j++){
	// 		for(k=0;k<Pdim;k++){
	// 			/* C(i,j) = sum(over k) A(i,k) * B(k,j) */
	// 			*(C+(i*Ndim+j)) += *(A+(i*Ndim+k)) *  *(B+(k*Pdim+j));
	// 		}
	// 	}
	// }

    // gather all local C matrixes across ranks in C matrix
    MPI_Gather(local_C, local_C_values, MPI_DOUBLE,
               C, local_C_values, MPI_DOUBLE, root, world_comm);
    
    //cleanup local A and C 
    free(local_A); free(local_C);


    /* ========================= END AREA OF INTEREST ========================== */
	run_time = MPI_Wtime() - start_time;
 
	if (rank == 0) {
        double max_time;
        MPI_Reduce(&run_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        printf("Order %d multiplication in %f seconds\n", order, run_time);

        dN = (double)order;
        mflops = 2.0 * dN * dN * dN / (1.0e6 * run_time);
        printf("Order %d multiplication at %f mflops\n", order, mflops);
        // Check the result
        cval = Pdim * AVAL * BVAL;
        errsq = 0.0;
        for (i = 0; i < Ndim; i++)
            for (j = 0; j < Mdim; j++) {
                err = C[i * Mdim + j] - cval;
                errsq += err * err;
            }

        if (errsq > TOL)
            printf("Errors in multiplication: %f\n", errsq);
        else
            printf("Hey, it worked!\n");

        free(A);
        free(C);
    } else {
        MPI_Reduce(&run_time, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    free(B);

    MPI_Finalize();
}