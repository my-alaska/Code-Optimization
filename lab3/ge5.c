//requires additional changes to the code to make it work

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <x86intrin.h>

#define max(a,b) ((a>b)?a:b)
#define BLCKSZ 8
#define _A_(i,j) A[i*SIZE+j]
#define IDX(i,j,size) (i*size+j)

static double gtod_ref_time_sec = 0.0;

/* Adapted from the bl2_clock() routine in the BLIS library */

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

int ge(double * A, int SIZE)
{
    register int i,j,k;
    register double division, Akk;
    for (k = 0; k < SIZE; k++) {
        Akk = _A_(k,k);
        for (i = k+1; i < SIZE; i++) {
            division = _A_(i,k)/Akk;
            for(j=k+1;j<SIZE;){
                if (j<max(SIZE-BLCKSZ, 0)){
                    _A_(i,j) = _A_(i,j)- _A_(k,j) *division;
                    _A_(i,j+1) = _A_(i,j+1)- _A_(k,j+1)*division;
                    _A_(i,j+2) = _A_(i,j+2)- _A_(k,j+2)*division;
                    _A_(i,j+3) = _A_(i,j+3)- _A_(k,j+3)*division;
                    _A_(i,j+4) = _A_(i,j+4)- _A_(k,j+4)*division;
                    _A_(i,j+5) = _A_(i,j+5)- _A_(k,j+5)*division;
                    _A_(i,j+6) = _A_(i,j+6)- _A_(k,j+6)*division;
                    _A_(i,j+7) = _A_(i,j+7)- _A_(k,j+7)*division;
                    j += BLCKSZ;
                }else{
                    _A_(i,j) = _A_(i,j)-_A_(k,j)*division;
                    j++;
                }
            }
        }
    }
    return 0;
}

int main( int argc, const char* argv[] )
{
    register int i,j,k;
    int iret;
    double dtime;
    int SIZE = 1500;

    double *M = malloc(SIZE*SIZE*sizeof(double));
    double **matrix = malloc(SIZE*sizeof(double*));
    for(int idx = 0; idx < SIZE; idx++) matrix[idx] = &(M[SIZE*idx]);

//  double matrix[SIZE][SIZE]; // TODO - make near optimal dynamic allocation

    srand(1);
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            matrix[i][j] = rand();
        }
    }
    printf("call GE");
    dtime = dclock();
    iret = ge(M, SIZE);
    dtime = dclock()-dtime;
    printf("Time: %le \n", dtime);

    double check=0.0;
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            check = check + matrix[i][j];
        }
    }
    printf( "Check: %le \n", check);
    fflush( stdout );


    return iret;
}
