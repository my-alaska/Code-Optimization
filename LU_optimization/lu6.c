#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define max(a,b) ((a>b)?a:b)
#define BLKSZ 8

// define macro for array access
#define IDX(i,j) (i*SIZE+j)

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

int LUPDecompose(double *A, register int SIZE, double Tol, int *P) {
    register int i, j, k, imax;
    register double maxA, ptr, absA;
    register double divisor, multiplier;
    int psize = SIZE;

    for (i = 0; i < SIZE; i++) P[i] = i;

    for (i = 0; i < SIZE; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < SIZE; k++)
            if ((absA = fabs(A[IDX(k,i)])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0;

        if (imax != i) {
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            // At some point we'll have to rewrite the matrix anyway so we might as well try to do it here
            // (It can be solved differently as shown in the next example)
            for(j = 0; j < SIZE; j++){
                ptr = A[IDX(i,j)];
                A[IDX(i,j)] = A[IDX(imax,j)];
                A[IDX(imax,j)] = ptr;
            }

            psize++;
        }

        // Use macro indexing

        divisor = A[IDX(i,i)];
        for (j = i + 1; j < SIZE; j++) {
            multiplier = A[IDX(j,i)]/divisor;
            A[IDX(j,i)] = multiplier;

            for (k=i+1;k<SIZE;){
                if (k < max(SIZE - BLKSZ, 0)){
                    A[IDX(j,k)] -= multiplier * A[IDX(i,k)];
                    A[IDX(j,k+1)] -= multiplier * A[IDX(i,k+1)];
                    A[IDX(j,k+2)] -= multiplier * A[IDX(i,k+2)];
                    A[IDX(j,k+3)] -= multiplier * A[IDX(i,k+3)];
                    A[IDX(j,k+4)] -= multiplier * A[IDX(i,k+4)];
                    A[IDX(j,k+5)] -= multiplier * A[IDX(i,k+5)];
                    A[IDX(j,k+6)] -= multiplier * A[IDX(i,k+6)];
                    A[IDX(j,k+7)] -= multiplier * A[IDX(i,k+7)];
                    k += BLKSZ;
                }else{
                    A[IDX(j,k)] -= multiplier * A[IDX(i,k)];
                    k++;
                }
            }
        }
    }

    P[SIZE] = psize;
    return 1;
}

int main( int argc, const char* argv[] ){
    register int i; // Removed j
    int iret;
    double dtime;
    int SIZE = 1500;
    double tol = 0.0000001;

    double *matrix = malloc(SIZE*SIZE*sizeof(double));
    int *P = malloc((SIZE+1)*sizeof(int));

    srand(1);
    for (i = 0; i < SIZE*SIZE; i++)  matrix[i] = rand(); // Why use 2 loops when we can use one?
//    printf("call LU decomposition");
    dtime = dclock();
    iret = LUPDecompose(matrix, SIZE, tol, P);
    dtime = dclock()-dtime;
//    printf("Time: ");
    printf("%.6f ", dtime);

    double check=0.0;
    for (i = 0; i < SIZE*SIZE; i++) check += matrix[i]; // Why use 2 loops when we can use one? (yep again)
//    printf( "\nCheck: %le \n", check);
    fflush( stdout );

    free(P);
    free(matrix);

    return iret;
}