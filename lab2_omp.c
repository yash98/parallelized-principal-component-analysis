#include <malloc.h>
#include <omp.h>
#include <math.h>

static float D_HAT_ARR[100000000];
static float U_arr[100000000];
static float VT_arr[100000000];

static float check_D[100000000];
static float check_usig[100000000];

void printMat(float* mat, int m, int n) {
    // mat is mxn
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            printf("%f", (*(mat+i*n+j)));
            if (j==n-1) {
                printf("\n");
            } else {
                printf(", ");
            }
        }
    }
}

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
        // printf("q\n");
        // printMat(q, m, n);

        // subtract previous direction components to make orthogonal 
        for (int j=0; j<n; j++) {
            if (j<i) {
                float dotp = subtractProjOnUnitVecFromRes(q+i,q+j,d+i,m,n);
                *(r+j*n+i) = dotp;
            } else {
                *(r+j*n+i) = 0.0;
            }
            // printf("q\n");
            // printMat(q, m, n);
        }

        // normalize
        float norm = 0.0;
        for (int j=0; j<m; j++) {
            norm += (*(q+j*n+i))*(*(q+j*n+i));
        }
        norm = sqrt(norm);

        *(r+i*n+i) = norm;

        for (int j=0; j<m; j++) {
            *(q+j*n+i) /= norm;
        }
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

void transpose(float* a, float* b, int m, int n) {
    // a = bT, a is mxn
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            *(a+i*n+j) = *(b+j*m+i);
        }
    }
}

void SVD_internal(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    // float* DT_D = (float*) malloc(sizeof(float)*N*N);
    float DT_D_arr[N*N];
    float* DT_D = DT_D_arr;
    
    // MT*M
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
           *(DT_D+i*N+j) = 0;
            for (int k=0; k<M; k++) {
                *(DT_D+i*N+j) += (*(D+k*N+i))*(*(D+k*N+j));
            }
        }
    }
    printf("DT_D\n");
    printMat(DT_D, N, N);

    // Eigen value and vectors
    // float* qi = (float*) malloc(sizeof(float)*N*N);
    // float* ri = (float*) malloc(sizeof(float)*N*N);
    // // float* di_even =(float*) malloc(sizeof(float)*N*N);
    // float* di_even = DT_D;
    // float* di_odd = (float*) malloc(sizeof(float)*N*N);
    // float* ei_even = (float*) malloc(sizeof(float)*N*N);
    // float* ei_odd = (float*) malloc(sizeof(float)*N*N);

    float qi_arr[N*N];
    float ri_arr[N*N];
    // float di_even_arr[N*N];
    float di_odd_arr[N*N];
    float ei_even_arr[N*N];
    float ei_odd_arr[N*N];

    float* qi = qi_arr;
    float* ri = ri_arr;
    // float* di_even = di_even_arr;
    float* di_even = DT_D;
    float* di_odd = di_odd_arr;
    float* ei_even = ei_even_arr;
    float* ei_odd = ei_odd_arr;

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
    // float error_even = 1000000000.0;
    // float error_odd = 0.0;
    float error =0.0;

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
            // errorI = &error_even;
            // errorIplus1 = &error_odd;
            errorIplus1 = &error;
        } else {
            di = di_odd;
            diPlus1 = di_even;
            ei = ei_odd;
            eiPlus1 = ei_even;
            // errorI = &error_odd;
            // errorIplus1 = &error_even;
            errorIplus1 = &error;
        }

        // all necessary calculation
        QRfactors(di, qi, ri, N, N);
        matmul(ri, qi, diPlus1, N, N, N);
        matmul(ei, qi, eiPlus1, N, N, N);
        // printf("q%d\n", numIter);
        // printMat(qi, N, N);
        // printf("r%d\n", numIter);
        // printMat(ri, N, N);
        // printf("e%d\n", numIter+1);
        // printMat(eiPlus1, N, N);
        // printf("d%d\n", numIter+1);
        // printMat(diPlus1, N, N);

        // Di+1 - D error stops decreasing
        *errorIplus1 = 0.0;
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                // float diff = *(diPlus1+i*N+j) - *(DT_D+i*N+j);
                float diff = *(diPlus1+i*N+j) - *(di+i*N+j);
                *errorIplus1 += diff*diff;
                // *errorIplus1 += diff;
                diff = *(eiPlus1+i*N+j) - *(ei+i*N+j);
                *errorIplus1 += diff*diff;
            }
        }
        numIter++;
        // if (*errorI<=*errorIplus1) {
        if (*errorIplus1<=0.00001) {
            break;
        }
    }
    printf("numIteration:%d\n", numIter);
    printf("Dinf\n");
    printMat(diPlus1, N, N);
    printf("Einf\n");
    printMat(eiPlus1, N, N);

    // float* sigmaInv = (float*) malloc(sizeof(float)*N*N);
    // float* v = (float*) malloc(sizeof(float)*N*N);
    float sigmaArr[M*N];
    float* sigmaPtr = sigmaArr;

    float sigmaInv_arr[N*M];
    float v_arr[N*N];
    float* sigmaInv = sigmaInv_arr;
    float* v = v_arr;

    for (int i=0; i<N; i++) {
        int newIndex = 0;
        float current = fabs(*(diPlus1+i*N+i));

        // calculate index
        for (int j=0; j<N; j++) {
            float other = fabs(*(diPlus1+j*N+j));
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
                *(sigmaInv+newIndex*M+newIndex) = 1.0/sqrt(current);
                // *(*(SIGMA)+newIndex*N+newIndex) = current;
                *(sigmaPtr+newIndex*N+newIndex) = sqrt(current);
                *(*(SIGMA)+newIndex) = sqrt(current);
            } else {
                *(sigmaInv+j*M+newIndex) = 0.0;
                // *(*(SIGMA)+j*N+newIndex) = 0.0;
                *(sigmaPtr+j*N+newIndex) = 0.0;
            }
        }
    }

    for (int i=0; i<N; i++) {
        for (int j=N; j<M; j++) {
            *(sigmaInv+i*M+j) = 0.0;
            *(sigmaPtr+j*N+i) = 0.0;
        }
    }

    printf("V\n");
    printMat(v, N, N);
    printf("V_T\n");
    printMat(*V_T, N, N);
    printf("SIGMA\n");
    printMat(*SIGMA, 1, N);
    printf("sigma\n");
    printMat(sigmaPtr, M, N);
    printf("sigmaInv\n");
    printMat(sigmaInv, N, M);

    // float* mv = (float*) malloc(sizeof(float)*M*N);
    float mv_arr[M*N];
    float* mv = mv_arr;

    matmul(D, v, mv, M, N, N);
    matmul(mv, sigmaInv, *U, M, N, M);

    // printf("sigma\n");
    // printMat(*SIGMA, 1, N);
    printf("U\n");
    printMat(*U, M, M);

    // checking
    float* checkD = check_D;
    float* checkUSig = check_usig;

    matmul(*U, sigmaPtr, checkUSig, M, M, N);
    matmul(checkUSig, *V_T, checkD, M, N, N);
    printf("checkUSig\n");
    printMat(checkUSig, M, N);
    printf("checkD\n");
    printMat(checkD, M, N);

    float checkError = 0.0;
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            float diff = (*(D+i*N+j))-(*(checkD+i*N+j));
            checkError += diff*diff;
        }
    }
    printf("checkError=%f\n", checkError);

    // free(DT_D);
    // // free(di_even);
    // free(di_odd);
    // free(qi);
    // free(ri);
    // free(ei_odd);
    // free(ei_even);
    // free(sigmaInv);
    // free(v);
    // free(mv);
}

void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T) {
    float* U1 = U_arr;
    float* VT1 = VT_arr;
    SVD_internal(M, N, D, &U1, SIGMA, &VT1);
    // U = (VT1)T, V_T = (U1)T
    transpose(*U, VT1, N, N);
    transpose(*V_T, U1, M, M);
}

void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    float totalEigenSum = 0.0;
    for (int i=0; i<N; i++) {
        totalEigenSum += *(SIGMA+i);
    }

    float eigenSum = 0.0;
    for (*K=0; *K<N; (*K)++) {
        eigenSum += *(SIGMA+(*K));
        float percent = eigenSum/totalEigenSum;
        if (percent >= retention/100.0) {
            (*K)++;
            break;
        }
    }

    // float* D_HAT_TEMP = (float*) calloc(M*(*K), sizeof(float));
    // float* D_HAT_TEMP = (float*) malloc(sizeof(float)*M*(*K));
    float* D_HAT_TEMP = D_HAT_ARR;

    // D_HAT = D*W, W = U[:*K]
    for (int i=0; i<M; i++) {
        for (int j=0; j<(*K); j++) {
            *(D_HAT_TEMP+i*(*K)+j) = 0.0;
            for (int k=0; k<N; k++) {
                *(D_HAT_TEMP+i*(*K)+j) += (*(D+i*N+k))*(*(U+k*N+j));
            }
        }
    }

    *D_HAT = D_HAT_TEMP;
    printf("K=%d\n", *K);
    printf("D_HAT\n");
    printMat(*D_HAT, M, *K);
}
