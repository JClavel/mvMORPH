/*---------Matrix parameterization for shared eigenvectors comparison-----------------*/
/*---------through the use of spherical Cholesky parameterization---------------------*/
/*-mvMORPH 1.0.5 - 2015 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/

#include "mvmorph.h"



SEXP spherical(SEXP param, SEXP variance, SEXP dim){
    int i, j, p, ind, index, f, col_start;
    char transa = 'T', transb = 'N';
    double one = 1.0 , zero = 0.0;
    p = INTEGER(dim)[0];

    SEXP U = PROTECT(allocMatrix(REALSXP,p,p));
    SEXP R = PROTECT(allocMatrix(REALSXP,p,p));
    SEXP V = PROTECT(allocMatrix(REALSXP,p,p));
    PROTECT(coerceVector(param,REALSXP));
    PROTECT(coerceVector(variance,REALSXP));
    // define pointers
    double *x = REAL(param), *upt = REAL(U), *corrMat = REAL(R), *varMat = REAL(V), *var_val = REAL(variance);
    // before the loop we fix the index of the param list
    upt[0]=1; // Fixed to 1 for computing the cholesky factor of a correlation matrix
    ind=0;
    col_start=0;
    index=0;
    
    // Compute the Cholesky factor through spherical parameterization (See Pinheiro & Bates 1996)
    for(i=1; i<p; i++){
        for(j=0; j<=i; j++){
            upt[i*p+j]=1;
        
            if(i==j){
                for(f=col_start; f<=(col_start+ind); f++) upt[i*p+j]*=sin(x[f]);
            }
        
            if(i!=j){
                upt[i*p+j]*=cos(x[col_start+j]);
                for(f=col_start; f<(col_start+j); f++) upt[i*p+j]*=sin(x[f]);
                upt[j*p+i]=0;
                index++;
            }

        }
        col_start=index;
        ind++;
    }
    
    // Do the matrix crossproduct of the Cholesky factors U to compute the correlation matrix R
    F77_CALL(dgemm)(&transa,&transb,&p,&p,&p,&one,REAL(U),&p,REAL(U),&p,&zero,REAL(R),&p FCONE FCONE);
    
    // Scale the correlation matrix R with the variance terms to obtain V
    for(i=0;i<p;i++){
        for(j=0;j<p;j++){
            varMat[i*p+j] = var_val[i]*var_val[j]*corrMat[i*p+j];
        }
    }
    
    UNPROTECT(5);
    return V;
}









