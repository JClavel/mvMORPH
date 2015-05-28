/* ---------- Orthogonal matrix using Givens rotations ----------- */
/* Julien Clavel - clavel@biologie.ens.fr/julien.clavel@hotmail.fr */
/* mvMORPH 1.0.3 --2015------------------------------------------- */

// Compute an orthogonal matrix using Givens rotations
// Use for SVD parametrization using eigenvectors and eigenvalues
// Avoid trigonometric functions following Golub & Van Loan 2013 and Stewart 1976
#include "mvmorph.h"


// Get sine and cosine of the Givens angle without trigonometrics functions
// Stewart (1976) - 3 / Golub & Van Loan (2013) - 5.1.10
/*static void getZ(double *p, double *s, double *c){
    // Avoid undefined values while using the non-trigonometric function
    if(p[0]>0.5 & p[0]<2 & p[0]!=1){
        p[0]=p[0]/4;
    }
    if(p[0]==1.0){
        c[0]=0;
        s[0]=1;
    }else if(abs(p[0])<1.0){
        s[0]=2*p[0];
        c[0]=sqrt(1-s[0]*s[0]);
    }else{
        c[0]=2/p[0];
        s[0]=sqrt(1-c[0]*c[0]);
    }
}*/

// Get sine and cosine for the Givens angle
static void getZ(double *p, double *s, double *c){
    c[0]=cos(p[0]);
    s[0]=sin(p[0]);
}

// Golub & Van Loan (2013) p. 241
static void updateA(double *A, double *c, double *s, int *ci, int *ck, int *ndim){
    int j, i=*ci, k=*ck, n=*ndim;
    double t1, t2;
    
    for(j=0; j<n; j++){
        t1=A[i*n+j];
        t2=A[k*n+j];
        A[j+n*i]=c[0]*t1-s[0]*t2;
        A[j+n*k]=s[0]*t1+c[0]*t2;
    }
}


// Function to compute an orthogonal matrix with n*(n-1)/2 Givens rotations
SEXP givens_ortho (SEXP Q, SEXP angle, SEXP ndim) {
  int n, i, j, ind_p=0;
  n = INTEGER(ndim)[0];
  SEXP c = PROTECT(allocVector(REALSXP,1));
  SEXP s = PROTECT(allocVector(REALSXP,1));
  SEXP rho = PROTECT(allocVector(REALSXP,1));
  SEXP p = PROTECT(coerceVector(angle,REALSXP));
  SEXP A = PROTECT(isReal(Q) ? duplicate(Q): coerceVector(Q, REALSXP));
  // NB: Q is expected to be an identity matrix on entry
  //SEXP A = PROTECT(allocMatrix(REALSXP,n,n));
  // transform A as an identity matrix
  //  for(i=0; i<n; i++){
  //      A[i+i*n]=1;
  //  }
    
    
  // update A
    for(i=0; i<(n-1); i++){
        for(j=(i+1); j<n; j++){
            REAL(rho)[0]=REAL(p)[ind_p];
            // Givens rotations
            getZ(REAL(rho),REAL(s),REAL(c));
            updateA(REAL(A),REAL(c),REAL(s),&i,&j,&n);
            ind_p++;
        }
    }
    
  UNPROTECT(5);
  return A;
}
