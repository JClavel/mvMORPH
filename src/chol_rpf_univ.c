/* Solve linear system with RPF Cholesky - Julien Clavel - mvMORPH 1.0.3 - 2014 */
/* Fast computation of the Cholesky factorization using BLAS 3 and half memmory */
/* storage (Rectangular Packed Format). The Cholesky factor is the used to      */
/* solve the generalized least square system. "row-major order version"         */

#include "mvmorph.h"


/* Fonction pour calculer le log-determinant d'une matrice stockée au format RPF "column-major order" */

// function to compute the determinant of the matrix A
static void determinant(double *det, double *ARF, int *n){
    int i, mod, n1, n2, n3, nn, na=*n;
    
    // Calcul du determinant
    det[0]=0;
    // Evalue si la matrice est d'ordre pair ou impair
    mod = na%2;
    // taille de la matrice A, N est pair
    if(mod == 0){
        
        // debut
        nn=0;
        n1=na/2;
        n2=na+1;
        n3=n1+1;
        for(i = 0; i<n1; i++){
            det[0]+=log(ARF[i+n1+nn]);
            det[0]+=log(ARF[i+n3+nn]);
            nn+=n2;
        }
        det[0]=det[0]*2.0;
        // taille de la matrice A est impair
    }else{
        nn = 0;
        n1 = na/2;
        n2 = na-n1;
        for(i = 0; i<n2; i++){
            if(i<n1){
                det[0]+=log(ARF[i+n1+nn]);
                det[0]+=log(ARF[i+n2+nn]);
            }else{
                det[0]+=log(ARF[i+n1+nn]);
            }
            nn+=na;
        }
        det[0]=det[0]*2.0;
    }
}

// Add measurement error to the matrix diagonal
static void ms_error(double *A, double *mserr, int *n){
    int na = *n, i, j, nn, n1, n2, mod;
  
    mod=na%2;
    if(mod == 0){
        nn=0;
        n1=na/2;
        n2=na+1;
        for(i=n1; i<na; i++){
            A[i+nn]+=mserr[i];
            nn+=n2;
        }
        nn=n1+1;
        for(j=0; j<n1; j++){
            A[j+nn]+=mserr[j];
            nn+=n2;
        }
    }else{ //N is odd
        nn = 0;
        n1 = na/2;
        n2 = na-n1;
        for(i = n1; i<na; i++){
            A[i+nn]+=mserr[i];
            nn+=na;
        }
        nn=n2;
        for(j = 0; j<n1; j++){
            A[j+nn]+=mserr[j];
            nn+=na;
        }
    }
}


// Factorisation de Cholesky RPF
SEXP   Chol_RPF_univ(SEXP A, SEXP D, SEXP dat, SEXP nterm, SEXP ndimA, SEXP mserr, SEXP ismserr){
    int n, info = 0, nt, err, one = 1;
    double alpha = 1.;
    char up = 'U', trans = 'T', diag = 'N', side = 'L', norm = 'N';
    nt = INTEGER(nterm)[0];
    n = INTEGER(ndimA)[0];
    err = INTEGER(ismserr)[0];
    PROTECT(A = coerceVector(A,REALSXP));
    PROTECT(mserr = coerceVector(mserr,REALSXP));
    SEXP Ddat = PROTECT(isReal(dat) ? duplicate(dat): coerceVector(dat, REALSXP));
    SEXP DD = PROTECT(isReal(D) ? duplicate(D): coerceVector(D, REALSXP));
    SEXP det = PROTECT(allocVector(REALSXP,1));
    
    // add measurement error
    if(err==1){
        ms_error(REAL(A),REAL(mserr), &n);
    }

    // decomposition de Cholesky
    F77_CALL(dpftrf)(&norm,&up,&n, REAL(A),&info);
    if (info != 0) {
        if (info > 0) error("the leading minor of order %d is not positive definite",info);
        error("argument %d of Lapack routine %s had invalid value",-info, "dpftrf");
    }
    // systeme lineaire U'x=D
    F77_CALL(dtfsm)(&norm, &side, &up, &trans, &diag, &n, &nt, &alpha, REAL(A), REAL(DD), &n);
    // systeme lineaire U'x=dat
    F77_CALL(dtfsm)(&norm, &side, &up, &trans, &diag, &n, &one, &alpha, REAL(A), REAL(Ddat), &n);
    
    // Calcul du determinant
    determinant(REAL(det),REAL(A),&n);
    
    // Liste: determinant, cholesky, X
    SEXP vec = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(vec, 0, A);
    SET_VECTOR_ELT(vec, 1, det);
    SET_VECTOR_ELT(vec, 2, DD);
    SET_VECTOR_ELT(vec, 3, Ddat);
    
    UNPROTECT (6);
    return (vec); 
}


// Factorisation de Cholesky RPF
SEXP   Chol_RPF_univ_only(SEXP A, SEXP ndimA, SEXP mserr, SEXP ismserr){
    int n, info = 0, err;
    char up = 'U', norm = 'N';
    
    n = INTEGER(ndimA)[0];
    err = INTEGER(ismserr)[0];
    PROTECT(A = coerceVector(A,REALSXP));
    PROTECT(mserr = coerceVector(mserr,REALSXP));
    SEXP det = PROTECT(allocVector(REALSXP,1));
    
    // add measurement error
    if(err==1){
        ms_error(REAL(A),REAL(mserr), &n);
    }
    
    // decomposition de Cholesky
    F77_CALL(dpftrf)(&norm,&up,&n, REAL(A),&info);
    if (info != 0) {
        if (info > 0) error("the leading minor of order %d is not positive definite",info);
        error("argument %d of Lapack routine %s had invalid value",-info, "dpftrf");
    }
    
    // Calcul du determinant
    determinant(REAL(det),REAL(A),&n);
    
    // Liste: determinant, cholesky, X
    SEXP vec = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(vec, 0, A);
    SET_VECTOR_ELT(vec, 1, det);
    
    UNPROTECT (4);
    return (vec);
}


// function to compute the quadratic product
SEXP   Chol_RPF_quadprod_column(SEXP U, SEXP resid, SEXP nterm){
    int n, info = 0, one = 1;
    double alpha = 1.;
    char up = 'U', trans = 'T', diag = 'N', side = 'L', norm = 'N';
    n = INTEGER(nterm)[0];
    PROTECT(U = coerceVector(U,REALSXP));
    SEXP Ddat = PROTECT(isReal(resid) ? duplicate(resid): coerceVector(resid, REALSXP));
    SEXP Bet = PROTECT(allocVector(REALSXP,1));
    double *beta = REAL(Bet), *data = REAL(Ddat), *chol = REAL(U);
    // systeme lineaire U'x=dat
    F77_CALL(dtfsm)(&norm, &side, &up, &trans, &diag, &n, &one, &alpha, chol, data, &n);
    if (info != 0){
        error("the %d argument had an illegal value",info);
    }
    
    // initialize
    beta[0]=0;
    int i = 0, round = down(n,4);
    // loop unrolling
    for(; i<round; i+=4){
        beta[0]+= square(i);
        beta[0]+= square(i+1);
        beta[0]+= square(i+2);
        beta[0]+= square(i+3);
    }
    for(; i<n; i++){
        beta[0]+= square(i);
    }
    
    UNPROTECT (3);
    return (Bet); 
}