#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

// define macros for size of the block and max function
#define max(a,b) ((a>b)?a:b)
#define BLKSZ 8

static double gtod_ref_time_sec = 0.0;
double dclock(){
    double the_time, norm_sec;
    struct timeval tv;
    gettimeofday( &tv, NULL );
    if ( gtod_ref_time_sec == 0.0 )
        gtod_ref_time_sec = ( double ) tv.tv_sec;
    norm_sec = ( double ) tv.tv_sec - gtod_ref_time_sec;
    the_time = norm_sec + tv.tv_usec * 1.0e-6;
    return the_time;
}

int LUPDecompose(double **A, register int SIZE, double Tol, int *P) {
    register int i, j, k, imax;
    register double maxA, *ptr, absA;
    register double divisor, multiplier;
    int psize = SIZE;

    for (i = 0; i < SIZE; i++) P[i] = i;

    for (i = 0; i < SIZE; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < SIZE; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0;

        if (imax != i) {
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            psize++;
        }


        divisor = A[i][i];
        for (j = i + 1; j < SIZE; j++) {
            multiplier = A[j][i]/divisor;
            A[j][i] = multiplier;

            // unroll the deepest loop to BLKSZ iterations
            for (k=i+1;k<SIZE;){
                if (k < max(SIZE - BLKSZ, 0)){
                    A[j][ k ] -= multiplier * A[i][ k ];
                    A[j][k+1] -= multiplier * A[i][k+1];
                    A[j][k+2] -= multiplier * A[i][k+2];
                    A[j][k+3] -= multiplier * A[i][k+3];
                    A[j][k+4] -= multiplier * A[i][k+4];
                    A[j][k+5] -= multiplier * A[i][k+5];
                    A[j][k+6] -= multiplier * A[i][k+6];
                    A[j][k+7] -= multiplier * A[i][k+7];
                    k += BLKSZ;
                }else{
                    A[j][k] -= multiplier * A[i][k];
                    k++;
                }
            }

        }
    }

    P[SIZE] = psize;
    return 1;
}

int main( int argc, const char* argv[] ){
    register int i,j;
    int iret;
    double dtime;
    int SIZE = 1500;
    double tol = 0.0000001;

    double *M = malloc(SIZE*SIZE*sizeof(double));
    double **matrix = malloc(SIZE*sizeof(double*));
    for(i = 0; i < SIZE; i++) matrix[i] = M + SIZE * i;

    int *P = malloc((SIZE+1)*sizeof(int));

    srand(1);
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            matrix[i][j] = rand();
        }
    }
//    printf("call LU decomposition");
    dtime = dclock();
    iret = LUPDecompose(matrix, SIZE, tol, P);
    dtime = dclock()-dtime;
//    printf("Time: ");
    printf("%.6f ", dtime);

    double check=0.0;
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            check = check + matrix[i][j];
        }
    }
//    printf( "\nCheck: %le \n", check);
    fflush( stdout );

    free(P);
    free(matrix);
    free(M);


    return iret;
}
