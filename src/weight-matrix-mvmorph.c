/*-Matrice de poids W pour un processus Ornstein-Uhlenbeck multivarie-----------------*/
/*-mvMORPH 1.0.3 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
/*-scale the weight matrix so that the rows sum to 1----------------------------------*/
/*-generalization to non-ultrametric trees and non-symmetric A matrix-----------------*/

#include "covar.h"



static void multi_weight_matrix (int *nchar, int *neps, double *epochs, double *lambda, double *S, double *S1, double *y) {
  double *elt;
  double t;
  int n = *nchar, np = *neps;
  int i, j, k, r;
  elt = Calloc(np*n,double);
  for (i = 0; i < np; i++) {
    t = epochs[0]-epochs[i];
    for (j = 0; j < n; j++)
      elt[i+np*j] = exp(-lambda[j]*t);
  }
    
  for (i = 0; i < np-1; i++) {
    for (j = 0; j < n; j++)
      elt[i+np*j] -= elt[i+1+np*j];
  }
    
    // calcul la matrice exponentielle en multipliant les vecteurs propres
  for (i = 0; i < np; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < n; k++) {
	y[j+n*(k+n*i)] = 0;
	for (r = 0; r < n; r++) 
	  y[j+n*(k+n*i)] += S[j+n*r]*elt[i+np*r]*S1[r+n*k]; // replace S[k+n*r] by S1[r+n*k] because the inverse of the eigenvector matrix S is explicitely given in S1
      }
    }
  }
    
  Free(elt);
}

// row standardization
static void row_stand(double *W, int *nt, int *ndim, int *ncol){
    int i, j, k, larg=*ndim, n=*nt, nc=*ncol;
    double sum_val;
   // la matrice W est en column major order, voir pour un unrolling?
    for(i=0; i<n; i++){
        sum_val=0;
        for(j=0; j<larg; j++){
            sum_val+=W[j*nc+i];
        }
        for(k=0; k<larg; k++){
            W[k*nc+i]/=sum_val;
        }
    }

}

// Construction de la matrice de variables indicatrices / building the design matrix
SEXP mvmorph_weights (SEXP nterm, SEXP epochs, SEXP lambda, SEXP S, SEXP S1, SEXP beta, SEXP root) {
  int nprotect = 0;
  SEXP W;
  double *wp, *y, *bp;
  int nchar, nt, *nreg, totreg, np, xdim[2], ptr, thetaO;
  int i, j, k, n, q;
  nchar = GET_LENGTH(lambda);  
  nt = INTEGER(nterm)[0];
  thetaO = INTEGER(root)[0];
  nreg = Calloc(nchar,int);
  totreg = 0;
    
    for (i = 0; i < nchar; i++) {
        nreg[i] = INTEGER(GET_DIM(VECTOR_ELT(VECTOR_ELT(beta,0),i)))[1];
        totreg += nreg[i];
    }
    
        xdim[0] = nt*nchar; xdim[1] = totreg;
        PROTECT(W = makearray(2,xdim)); nprotect++;
        for (i = 0; i < nt; i++) {
            np = GET_LENGTH(VECTOR_ELT(epochs,i));
            y = Calloc(nchar*nchar*np,double);
            multi_weight_matrix(&nchar,&np,REAL(VECTOR_ELT(epochs,i)),REAL(lambda),REAL(S),REAL(S1),y);
            
            for (n = 0, ptr = 0; n < nchar; ptr += nt*nchar*nreg[n++]) {
                wp = &(REAL(W))[ptr];
                bp = REAL(VECTOR_ELT(VECTOR_ELT(beta,i),n));
                for (j = 0; j < nchar; j++) {
                    for (k = 0; k < nreg[n]; k++) {
                        wp[i+nt*(j+nchar*k)] = 0;
                        for (q = 0; q < np; q++) {
                            wp[i+nt*(j+nchar*k)] += y[n+nchar*(j+nchar*q)]*bp[q+np*k];
                        }
                    }
                }
            }
            
            
            Free(y);
        }
        Free(nreg);
    
    // Row standardize the matrix if the root is not estimated
    if(thetaO==1){
        int nsp;
        nsp=nt*nchar;
        row_stand(REAL(W),&nsp,&totreg,&nsp);
    }
    
  UNPROTECT(nprotect);
  return W;
}

