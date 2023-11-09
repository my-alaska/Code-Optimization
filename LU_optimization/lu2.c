#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

static double gtod_ref_time_sec = 0.0;
double dclock()
{
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
    register int i, j, k, imax; //storing frequently used variables in a register
    register double maxA, absA;
    double *ptr;

    for (i = 0; i <= SIZE; i++)
        P[i] = i;

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

            P[SIZE]++;
        }

        for (j = i + 1; j < SIZE; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < SIZE; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;
}

int main( int argc, const char* argv[] ){
    register int i, j, SIZE; //storing frequently used variables in a register
    int iret;
    double dtime;
    double tol = 0.0000001;
    SIZE = 1500;

    double **matrix = malloc(SIZE*sizeof(double*));
    for(i = 0; i < SIZE; i++) matrix[i] = malloc(SIZE*sizeof(double));
    int *P = malloc(SIZE*sizeof(int));

    srand(1);
    for (i = 0; i < SIZE; i++) for (j = 0; j < SIZE; j++) matrix[i][j] = rand();

//    printf("call LU decomposition");
    dtime = dclock();
    iret = LUPDecompose(matrix, SIZE, tol, P);
    dtime = dclock()-dtime;
//    printf("Time: ");
    printf("%.6f ", dtime);

    double check=0.0;
    for (i = 0; i < SIZE; i++) for (j = 0; j < SIZE; j++) check = check + matrix[i][j];

//    printf( "Check: %le \n", check);
    fflush( stdout );

    free(P);
    for(i = 0; i < SIZE; i++)free(matrix[i]);
    free(matrix);


    return iret;
}