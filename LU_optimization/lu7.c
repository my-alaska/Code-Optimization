#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define max(a,b) ((a>b)?a:b)
#define BLKSZ 8
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
    register int i, j, k, imax, iIdx, jIdx;
    register double maxA, absA;
    register double divisor, multiplier;
    int psize = SIZE;

    for (i = 0; i < SIZE; i++) P[i] = i;

    for (i = 0; i < SIZE; i++) {
        // we use indexing function (P array) to access the row o the matrix so we don't need to rewrite it
        iIdx = P[i];
        maxA = 0.0;
        imax = iIdx;

        for (k = i; k < SIZE; k++)
            // In order to not rewrite parts of the matrix we can use the P array as our indexing function
            if ((absA = fabs(A[IDX(P[k],i)])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0;

        if (imax != i) {
            P[i] = P[imax];
            P[imax] = iIdx;
            iIdx = P[i];
            psize++;
        }

        divisor = A[IDX(iIdx,i)];
        for (j = i + 1; j < SIZE; j++) {
            // So as not to call the indexing function(accessing the array)
            // many time inside the loop we can do it once and store it in a registered index (iIdx and jIdx)
            jIdx = P[j];
            multiplier = A[IDX(jIdx,i)]/divisor;
            A[IDX(jIdx,i)] = multiplier;

            for (k=i+1;k<SIZE;){
                if (k < max(SIZE - BLKSZ, 0)){
                    A[IDX(jIdx,k)] -= multiplier * A[IDX(iIdx,k)];
                    A[IDX(jIdx,k+1)] -= multiplier * A[IDX(iIdx,k+1)];
                    A[IDX(jIdx,k+2)] -= multiplier * A[IDX(iIdx,k+2)];
                    A[IDX(jIdx,k+3)] -= multiplier * A[IDX(iIdx,k+3)];
                    A[IDX(jIdx,k+4)] -= multiplier * A[IDX(iIdx,k+4)];
                    A[IDX(jIdx,k+5)] -= multiplier * A[IDX(iIdx,k+5)];
                    A[IDX(jIdx,k+6)] -= multiplier * A[IDX(iIdx,k+6)];
                    A[IDX(jIdx,k+7)] -= multiplier * A[IDX(iIdx,k+7)];
                    k += BLKSZ;
                }else{
                    A[IDX(jIdx,k)] -= multiplier * A[IDX(iIdx,k)];
                    k++;
                }
            }
        }
    }

    P[SIZE] = psize;
    return 1;
}

int main( int argc, const char* argv[] ){
    register int i;
    int iret;
    double dtime;
    int SIZE = 1500;
    double tol = 0.0000001;

    double *matrix = malloc(SIZE*SIZE*sizeof(double));
    int *P = malloc((SIZE+1)*sizeof(int));

    srand(1);
    for (i = 0; i < SIZE*SIZE; i++)  matrix[i] = rand();
//    printf("call LU decomposition");
    dtime = dclock();
    iret = LUPDecompose(matrix, SIZE, tol, P);
    dtime = dclock()-dtime;
//    printf("Time: ");
    printf("%.6f ", dtime);

    double check=0.0;
    for (i = 0; i < SIZE*SIZE; i++) check += matrix[i];
//    printf( "\nCheck: %le \n", check);
    fflush( stdout );

    free(P);
    free(matrix);

    return iret;
}