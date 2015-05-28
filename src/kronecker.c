/*-----------------------Kronecker product computation--------------------------------*/
/*-mvMORPH 1.0.3 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/

#include "mvmorph.h"


static void copy_sparse(int *IA, int *JA, int *nrow, double *restrict V, double *restrict A){
    
    int n, i, j, init, end, inc=0;
    n = *nrow;
#define V(r, c) (V[(r)*n + (c)])
    
    
    for( j = 0 ; j<n ; j++){
        init=IA[j];
        end=IA[j+1];
        for ( i = init; i < end; i++){
            A[inc]=V(j,JA[i]);
            inc++;
        }  
    }
    
    
}

/*
 static void sum_sparse(int *IA, int *JA, int *nrow, double *restrict V, double *restrict A){
    
    int n, i, j, init, end, inc=0;
    n = *nrow;
#define V(r, c) (V[(r)*n + (c)])
    
    for( j = 0 ; j<n ; j++){
        init=IA[j];
        end=IA[j+1];
        for ( i = init; i < end; i++){
            A[inc]+=V(j,JA[i]);
            inc++;
        }
    }
}
 */

// Adapted from Yadav & Michalak 2013
static void kronecker_eb(int *rrows, int *crows, double *R, double *C, double *V, double *beta){
    int i, j, k, l, m, n=0, x=0, rdim, cdim;
    double ebtrans;
    rdim=*rrows;
    cdim=*crows;
    
    for (i = 0; i <= rdim-1; i++){
        for ( j = (i*cdim); j <= (i*cdim)+(cdim-1); j++){
            m=i*rdim;
            for ( k = (j*rdim); k <= (j*rdim)+(rdim-1); k++){
                n=(j*cdim)-(i*x);
                for (l = (k*cdim); l <= (k*cdim)+(cdim-1); l++){
                    
                    if(beta[m]==0){
                        ebtrans=C[n];
                    }else{
                        ebtrans=(exp(beta[m]*C[n])-1)/beta[m];
                    }
                    
                    V[l] = R[m]*ebtrans;
                    n++;
                }
                m++;
            }
        }
        x=n;
    }
}

static void kronecker_eb_sum(int *rrows, int *crows, double *R, double *C, double *V, double *beta){
    int i, j, k, l, m, n=0, x=0, rdim, cdim;
    double ebtrans;
    rdim=*rrows;
    cdim=*crows;
    
    for (i = 0; i <= rdim-1; i++){
        for ( j = (i*cdim); j <= (i*cdim)+(cdim-1); j++){
            m=i*rdim;
            for ( k = (j*rdim); k <= (j*rdim)+(rdim-1); k++){
                n=(j*cdim)-(i*x);
                for (l = (k*cdim); l <= (k*cdim)+(cdim-1); l++){
                    
                    if(beta[m]==0){
                        ebtrans=C[n];
                    }else{
                        ebtrans=(exp(beta[m]*C[n])-1)/beta[m];
                    }
                    
                    V[l] += R[m]*ebtrans;
                    n++;
                }
                m++;
            }
        }
        x=n;
    }
}

static void kronecker_sum_eb_bm(int *rrows, int *crows, double *R1, double *R2, double *C1, double *C2, double *V, double *beta){
    int i, j, k, l, m, n=0, x=0, rdim, cdim;
    double ebtrans, tmp1, tmp2;
    rdim=*rrows;
    cdim=*crows;
    
    for (i = 0; i <= rdim-1; i++){
        for ( j = (i*cdim); j <= (i*cdim)+(cdim-1); j++){
            m=i*rdim;
            for ( k = (j*rdim); k <= (j*rdim)+(rdim-1); k++){
                n=(j*cdim)-(i*x);
                for (l = (k*cdim); l <= (k*cdim)+(cdim-1); l++){
                    
                    if(beta[m]==0){
                        ebtrans=C2[n];
                    }else{
                        ebtrans=(exp(beta[m]*C2[n])-1)/beta[m];
                    }
                    tmp1=R1[m]*ebtrans;
                    tmp2=R2[m]*C1[n];
                    
                    V[l] = tmp1+tmp2 ;
                    n++;
                }
                m++;
            }
        }
        x=n;
    }
}

static void kronecker_sum(int *rrows, int *crows, double *restrict R, double *restrict C, double *restrict V){
    int i, j, k, l, m, n=0, x=0, rdim, cdim;
    rdim=*rrows;
    cdim=*crows;
    
    for (i = 0; i <= rdim-1; i++){
        for ( j = (i*cdim); j <= (i*cdim)+(cdim-1); j++){
            m=i*rdim;
            for ( k = (j*rdim); k <= (j*rdim)+(rdim-1); k++){
                n=(j*cdim)-(i*x);
                for (l = (k*cdim); l <= (k*cdim)+(cdim-1); l++){
                    V[l] += R[m]*C[n];
                    n++;
                }
                m++;
            }
        }
        x=n;
    }
    
    
}

static void kronecker_C(int *rrows, int *crows, double *restrict R, double *restrict C, double *restrict V){
    int i, j, k, l, m, n=0, x=0, rdim, cdim;
    rdim=*rrows;
    cdim=*crows;
    
    for (i = 0; i <= rdim-1; i++){
        for ( j = (i*cdim); j <= (i*cdim)+(cdim-1); j++){
            m=i*rdim;
            for ( k = (j*rdim); k <= (j*rdim)+(rdim-1); k++){
                n=(j*cdim)-(i*x);
                for (l = (k*cdim); l <= (k*cdim)+(cdim-1); l++){
                    V[l] = R[m]*C[n];
                    n++;
                }
                m++;
            }
        }
        x=n;
    }
    
    
}

// for shift model; i.e. only two mapped matrices
// kronecker sum of BM and EB models
SEXP kronecker_shiftEB_BM(SEXP R1, SEXP R2, SEXP C1, SEXP C2, SEXP beta, SEXP Rrows, SEXP Crows){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R1 = coerceVector(R1,REALSXP)); nprotect++;
    PROTECT(R2 = coerceVector(R2,REALSXP)); nprotect++;
    PROTECT(C1 = coerceVector(C1,REALSXP)); nprotect++;
    PROTECT(C2 = coerceVector(C2,REALSXP)); nprotect++;
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    
    kronecker_sum_eb_bm(&rrows, &crows, REAL(R1), REAL(R2), REAL(C1), REAL(C2), REAL(V), REAL(beta));
    
    UNPROTECT(nprotect);
    return V;
}

// for shift model; i.e. only two mapped matrices
// generalized for EB and OU
SEXP kronecker_shiftEB_OU(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP V){
    int rrows, crows, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(V = coerceVector(V,REALSXP)); nprotect++;
    
    kronecker_eb_sum(&rrows, &crows, REAL(R), REAL(C), REAL(V), REAL(beta));
    
    UNPROTECT(nprotect);
    return V;
}


// Standard kronecker product + sum of dense VCV matrices
SEXP kronecker_shift(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP V){
    int rrows, crows, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(V = coerceVector(V,REALSXP)); nprotect++;
    
    kronecker_sum(&rrows, &crows, REAL(R), REAL(C), REAL(V));
    
    UNPROTECT(nprotect);
    return V;
}

// Standard Kronecker + Sum for sparse matrix
// shift OU - EB
SEXP kroneckerSpar_shift_EB_BM(SEXP R1, SEXP R2, SEXP C1, SEXP C2, SEXP beta, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R1 = coerceVector(R1,REALSXP)); nprotect++;
    PROTECT(R2 = coerceVector(R2,REALSXP)); nprotect++;
    PROTECT(C1 = coerceVector(C1,REALSXP)); nprotect++;
    PROTECT(C2 = coerceVector(C2,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;

    
    kronecker_sum_eb_bm(&rrows, &crows, REAL(R1), REAL(R2), REAL(C1), REAL(C2), REAL(V), REAL(beta));
    
    copy_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    
    UNPROTECT(nprotect);
    return A;
}

// Standard Kronecker + Sum for sparse matrix
// shift OU - EB
/*SEXP kroneckerSpar_shift_OU_EB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    
    kronecker_eb(&rrows, &crows, REAL(R), REAL(C), REAL(V), REAL(beta));
    
    sum_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    
    UNPROTECT(nprotect);
    return A;
}*/
// Sparse OU-EB
SEXP kroneckerSpar_shift_OU_EB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP V, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;

    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(V = coerceVector(V,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    
    kronecker_eb_sum(&rrows, &crows, REAL(R), REAL(C), REAL(V), REAL(beta));
    
    copy_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    
    UNPROTECT(nprotect);
    return A;
}

// Standard Kronecker + Sum for sparse matrix
/*SEXP kroneckerSpar_shift(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    
    kronecker_C(&rrows, &crows, REAL(R), REAL(C), REAL(V));
    
    sum_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    
    UNPROTECT(nprotect);
    return A;
    
}*/
// Sparse OU-BM
SEXP kroneckerSpar_shift(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP V, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(V = coerceVector(V,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    
    kronecker_sum(&rrows, &crows, REAL(R), REAL(C), REAL(V));
    
    copy_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    
    UNPROTECT(nprotect);
    return A;
}


// Standard kronecker product of dense VCV matrices
SEXP kronecker_mvmorph(SEXP R, SEXP C, SEXP Rrows, SEXP Crows){
    int rrows, crows, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
        
    kronecker_C(&rrows, &crows, REAL(R), REAL(C), REAL(V));
    
    UNPROTECT(nprotect);
    return V;
}


// kronecker product of VCV matrices in sparse format (need to be optimized)
SEXP kroneckerSumSpar(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP dimlist, SEXP IA, SEXP JA, SEXP A){
    int f, rrows, crows, ndim, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    ndim=INTEGER(dimlist)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,VECSXP)); nprotect++;
    PROTECT(C = coerceVector(C,VECSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
 
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    // Note that zeroing the vector is necessary because of the sum of list elements
    memset(REAL(V),0,(vdim*vdim)*sizeof(double));

    
    for(f=0; f<ndim; f++){// traverse la liste de vcv
        
        kronecker_sum(&rrows, &crows, REAL(VECTOR_ELT(R,f)), REAL(VECTOR_ELT(C,f)), REAL(V));
    
    }// End
    
        copy_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
        
        UNPROTECT(nprotect);
        return A;
  
}

// kronecker product of dense VCV matrices (need to be optimized)

SEXP kroneckerSum(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP dimlist){
    int f, rrows, crows, ndim, vdim, nprotect=0;
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    ndim=INTEGER(dimlist)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,VECSXP)); nprotect++;
    PROTECT(C = coerceVector(C,VECSXP)); nprotect++;
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    // Note that zeroing the vector is necessary because of the sum of list elements
    memset(REAL(V),0,(vdim*vdim)*sizeof(double));
    
    for(f=0; f<ndim; f++){// traverse la liste de vcv
        
        kronecker_sum(&rrows, &crows, REAL(VECTOR_ELT(R,f)), REAL(VECTOR_ELT(C,f)), REAL(V));
        
    }

    UNPROTECT(nprotect);
    return V;
}


// Kronecker product of matrices with ACDC transformation in Sparse format

SEXP kroneckerSparEB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A){
    int rrows, crows, vdim, nprotect=0;
  
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
    PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
    PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;

    kronecker_eb(&rrows, &crows, REAL(R), REAL(C), REAL(V), REAL(beta));
    copy_sparse(INTEGER(IA),INTEGER(JA),&vdim,REAL(V),REAL(A));
    // saving in sparse format need to be optimized
        UNPROTECT(nprotect);
        return A;
    }


// Kronecker product of matrices with ACDC transformation in Dense format

SEXP kroneckerEB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows){
    int rrows, crows, vdim, nprotect=0;
    
    rrows=INTEGER(Rrows)[0];
    crows=INTEGER(Crows)[0];
    vdim=rrows*crows;
    
    PROTECT(R = coerceVector(R,REALSXP)); nprotect++;
    PROTECT(C = coerceVector(C,REALSXP)); nprotect++;
    
    SEXP V = PROTECT(allocMatrix(REALSXP,vdim,vdim)); nprotect++;
    
    kronecker_eb(&rrows, &crows, REAL(R), REAL(C), REAL(V), REAL(beta));
 
    
    UNPROTECT(nprotect);
    return V;
}



