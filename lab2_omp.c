#include <malloc.h>
#include <omp.h>
#include <math.h>

float subtractProjOnUnitVecFromRes(float* result, float* unitVec, float* generalVec, int m, int n) {
    float dotProdRes = 0.0;
    // dotproduct
    for (int i=0; i<m; i++) {
        dotProdRes += (*(unitVec+i*n))*(*(generalVec+i*n));
    }

    // subtract component
    for (int i=0; i<m; i++) {
        *(result+i*n) -= dotProdRes*(*(unitVec+i*n));
    }
    return dotProdRes;
}

void QRfactors(float* d, float* q, float* r, int m, int n) {
    for (int i=0; i<n; i++) {
        // copy ai to ei
        for (int j=0; j<m; j++) {
            *(q+j*n+i) = *(d+j*n+i);
        }

        // subtract previous direction components to make orthogonal 
        for (int j=0; j<n; j++) {
            if (j<i) {
                float dotp = subtractProjOnUnitVecFromRes(q+i,q+j,d+i,m,n);
                *(r+j*n+i) = dotp;
            } else {
                *(r+j*n+i) = 0.0;
            }
        }

        // normalize
        float norm = 0.0;
        for (int j=0; j<m; j++) {
            norm += (*(q+j*n+i))*(*(q+j*n+i));
        }
        norm = sqrt(norm);
        *(r+i*n+i) = norm;
    }
}

void matmul(float* a, float* b, float*c, int l, int m, int n) {
    // c = a*b, a is lxm, b is mxn, cis lxn
    for (int i=0; i<l; i++) {
        for (int j=0; j<n; j++) {
            *(c+i*n+j) = 0.0;
            for (int k=0; k<m; k++) {
                *(c+i*n+j) += (*(a+i*m+k))*(*(b+k*n+j));
            }
        }
    }
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    float* DT_D = (float*) malloc(sizeof(float)*N*N);
    
    // MT*M
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
           *(DT_D+i*N+j) = 0;
            for (int k=0; k<M; k++) {
                *(DT_D+i*N+j) += (*(D+k*N+i))*(*(D+k*N+j));
            }
        }
    }

    // Eigen value and vectors
    float* qi = (float*) malloc(sizeof(float)*N*N);
    float* ri = (float*) malloc(sizeof(float)*N*N);
    // float* di_even =(float*) malloc(sizeof(float)*N*N);
    float* di_even = DT_D;
    float* di_odd = (float*) malloc(sizeof(float)*N*N);
    float* ei_even = (float*) malloc(sizeof(float)*N*N);
    float* ei_odd = (float*) malloc(sizeof(float)*N*N);

    // init D0 & E0
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (i==j) {
                *(ei_even+i*N+j) = 1.0;
            } else {
                *(ei_even+i*N+j) = 0.0;
            }
            // *(di_even+i*N+j) = *(DT_D+i*N+j);
        }
    }

    int numIter = 0;
    float error_even = 1000000000.0;
    float error_odd = 0.0;

    float* di;
    float* diPlus1;
    float* ei;
    float* eiPlus1;
    float* errorI;
    float* errorIplus1;
    // until convergence
    while (numIter<10000) {
        // set saved value for window 2
        if (numIter%2==0) {
            di = di_even;
            diPlus1 = di_odd;
            ei = ei_even;
            eiPlus1 = ei_odd;
            errorI = &error_even;
            errorIplus1 = &error_odd;
        } else {
            di = di_odd;
            diPlus1 = di_even;
            ei = ei_odd;
            eiPlus1 = ei_even;
            errorI = &error_odd;
            errorIplus1 = &error_even;
        }

        // all necessary calculation
        QRfactors(di, qi, ri, N, N);
        matmul(ri, qi, diPlus1, N, N, N);
        matmul(ei, qi, eiPlus1, N, N, N);

        // Di+1 - D error stops decreasing
        *errorIplus1 = 0.0;
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                // float diff = *(diPlus1+i*N+j) - *(DT_D+i*N+j);
                float diff = *(diPlus1+i*N+j) - *(di+i*N+j);
                *errorIplus1 += diff*diff;
            }
        }
        if (errorI<=errorIplus1) {
            break;
        }
        numIter++;
    }

    float* sigmaInv = (float*) malloc(sizeof(float)*N*N);
    float* v = (float*) malloc(sizeof(float)*N*N);

    for (int i=0; i<N; i++) {
        int newIndex = 0;
        float current = *(diPlus1+i*N+i);

        // calculate index
        for (int j=0; j<N; j++) {
            float other = *(diPlus1+j*N+j);
            if (other>current) {
                newIndex++;
            } else if (other==current) {
                if (i<j) {
                    newIndex++;
                }
            }
        }

        // paste according to index
        for (int j=0; j<N; j++) {
            *(*(V_T)+newIndex*N+j) = *(eiPlus1+j*N+newIndex);
            *(v+j*N+newIndex) = *(eiPlus1+j*N+newIndex);
            if (j==newIndex) {
                *(sigmaInv+newIndex*N+newIndex) = 1.0/current;
                *(*(SIGMA)+newIndex*N+newIndex) = current;
            } else {
                *(sigmaInv+j*N+newIndex) = 0.0;
                *(*(SIGMA)+j*N+newIndex) = 0.0;
            }
        }
    }

    float* mv = (float*) malloc(sizeof(float)*M*N);

    matmul(D, v, mv, M, N, N);
    matmul(mv, sigmaInv, *U, M, N, N);

    free(DT_D);
    // free(di_even);
    free(di_odd);
    free(qi);
    free(ri);
    free(ei_odd);
    free(ei_even);
    free(sigmaInv);
    free(v);
    free(mv);
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
