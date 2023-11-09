#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include <immintrin.h>

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
    register __m256d mm_multiplier;
    register __m256d tmp0, tmp1, tmp2, tmp3;
    int psize = SIZE;

    for (i = 0; i < SIZE; i++) P[i] = i;

    for (i = 0; i < SIZE; i++) {
        iIdx = P[i];
        maxA = 0.0;
        imax = iIdx;

        for (k = i; k < SIZE; k++)
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
            jIdx = P[j];
            multiplier = A[IDX(jIdx,i)]/divisor;
            A[IDX(jIdx,i)] = multiplier;

            mm_multiplier[3] = mm_multiplier[2] =  mm_multiplier[1] = mm_multiplier[0] = multiplier;

            // use 256 bit vector operation
            for (k=i+1;k<SIZE;){
                if (k < max(SIZE - BLKSZ, 0)){
                    tmp0 = _mm256_loadu_pd(A + IDX(jIdx,k));
                    tmp2 = _mm256_loadu_pd(A + IDX(jIdx,k+4));

                    tmp1 = _mm256_loadu_pd(A + IDX(iIdx,k));
                    tmp3 = _mm256_loadu_pd(A + IDX(iIdx,k+4));

                    tmp1 = _mm256_mul_pd(tmp1,mm_multiplier);
                    tmp3 = _mm256_mul_pd(tmp3,mm_multiplier);

                    tmp0 = _mm256_sub_pd(tmp0,tmp1);
                    tmp2 = _mm256_sub_pd(tmp2,tmp3);

                    _mm256_storeu_pd(A + IDX(jIdx,k),tmp0);
                    _mm256_storeu_pd(A + IDX(jIdx,k+4),tmp2);
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