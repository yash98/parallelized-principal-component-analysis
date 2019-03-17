#include <malloc.h>
#include <omp.h>
#include <math.h>

float subtractProjOnUnitVecFromRes(float* result, float* unitVec, float* generalVec, int m, int n) {
    float dotProdRes = 0.0;
    for (int i=0; i<m; i++) {
        dotProdRes += (*(unitVec+i*n))*(*(generalVec+i*n));
    }
    for (int i=0; i<m; i++) {
        *(result+i*n) -= dotProdRes*(*(unitVec+i*n));
    }
    return dotProdRes;
}

void QRfactors(float* d, float* q, float* r, int m, int n) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            *(q+j*n+i) = *(d+j*n+i);
        }
        for (int j=0; j<n; j++) {
            if (j<i) {
                float dotp = subtractProjOnUnitVecFromRes(q+i,q+j,d+i,m,n);
                *(r+j*n+i) = dotp;
            } else {
                *(r+j*n+i) = 0.0;
            }
        }

        float norm = 0.0;
        for (int j=0; j<m; j++) {
            norm += (*(q+j*n+i))*(*(q+j*n+i));
        }
        norm = sqrt(norm);
        *(r+i*n+i) = norm;
    }
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    float* DT_D = (float*) malloc(sizeof(float)*N*N);
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
           *(DT_D+i*N+j) = 0;
            for (int k=0; k<M; k++) {
                *(DT_D+i*N+j) += (*(D+k*N+i))*(*(D+k*N+j));
            }
        }
    }

    float* qi = (float*) malloc(sizeof(float)*N*N);
    float* ri = (float*) malloc(sizeof(float)*N*N);
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
