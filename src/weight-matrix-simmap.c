/*-Matrice de poids W pour un processus Ornstein-Uhlenbeck multivarie-----------------*/
/*-mvMORPH 1.0.2 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
/*-Modified with permission from A. King OUCH ----------------------------------------*/

#include "covar.h"

static void simmap_weight_matrix (int *nchar, int *neps, double *epochs, double *lambda, double *S, double *y) {
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
  for (i = 0; i < np; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < n; k++) {
	y[j+n*(k+n*i)] = 0;
	for (r = 0; r < n; r++) 
	  y[j+n*(k+n*i)] += S[j+n*r]*elt[i+np*r]*S[k+n*r];
      }
    }
  }
  Free(elt);
}

SEXP simmap_weights (SEXP nterm, SEXP epochs, SEXP lambda, SEXP S, SEXP beta) {
  int nprotect = 0;
  SEXP W;
  double *wp, *y, *bp;
  int nchar, nt, *nreg, totreg, np, xdim[2], ptr;
  int i, j, k, n, q;
  nchar = GET_LENGTH(lambda);  
  nt = INTEGER(nterm)[0];
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
    simmap_weight_matrix(&nchar,&np,REAL(VECTOR_ELT(epochs,i)),REAL(lambda),REAL(S),y);
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
  UNPROTECT(nprotect);
  return W;
}

