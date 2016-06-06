/*-Matrice de covariance pour un processus Ornstein-Uhlenbeck multivarie-version-OUCH-*/
/*-mvMORPH 1.0.2 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
/*-modified version of the covar-matrix.c code from OUCH package----------------------*/
#include "covar.h"


// stationary covariance for Complex eigenvalues
static void simmap_covar_matrix_complex (int *nchar,
                                         double *bt,
                                         Rcomplex *lambda_val,
                                         Rcomplex *S_val,
                                         Rcomplex *S1_val,
                                         double *sigmasq,
                                         int *nterm,
                                         double *V) {
    

    // complex version
    double complex *eltj, *elti, *W, *U, *tmp1, *lambda, *S, *S1;
    double sij, ti, tj, tmp2;
    int n = *nchar, nt = *nterm;
    int i, j, k, l, r, s;
    
    // alloc complex vectors
    U = Calloc(n*n,double complex);
    W = Calloc(n*n,double complex);
    tmp1 = Calloc(n*n,double complex);
    eltj = Calloc(n,double complex);
    elti = Calloc(n,double complex);
    S = Calloc(n*n,double complex);
    S1 = Calloc(n*n,double complex);
    lambda = Calloc(n,double complex);
    
    //zeroing vectors & transform to C complex structure
    for(i = 0; i<n; i++){
        lambda[i]=comp(lambda_val[i]);
        for(j =0; j<n; j++){
            S[i+j*n]=comp(S_val[i+j*n]);
            S1[i+j*n]=comp(S1_val[i+j*n]);
            U[i+j*n] = 0;
            W[i+j*n] = 0;
        }
    }
    
    // computing the P^-1%*%Sigma%*%P^-T
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    U[i+j*n] += S1[i+k*n]*sigmasq[k+l*n]*S1[j+l*n];
                }
            }
        }
    }
    
    // fill in the covariance
    for (i = 0; i < nt; i++) {
        for (j = 0; j <= i; j++) {
            ti = bt[i+i*nt];
            sij = bt[i+j*nt];
            tj = bt[j+j*nt];
            
            // complex exponential with time
            for (k = 0; k < n; k++) {
                elti[k] = cexp(-lambda[k]*(ti-sij));
                eltj[k] = cexp(-lambda[k]*(tj-sij));
            }
            
            // Integral parts
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    W[k+l*n] = elti[k]*U[k+l*n]*eltj[l]/(lambda[k]+lambda[l]);
                }
            }
            
            // computing the P%*%Sigma%*%P^T
            for (r = 0; r < n; r++) {
                for (s = 0; s < n; s++) {
                    tmp1[r+s*n] = 0;
                    for (k = 0; k < n; k++) {
                        for (l = 0; l < n; l++) {
                            tmp1[r+s*n] += S[r+k*n]*W[k+l*n]*S[s+l*n];
                        }
                    }
                }
            }

            
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    // Save the real parts
                    tmp2 = creal(tmp1[k+l*n]);
                    
                    V[i+nt*(k+n*(j+nt*l))] = tmp2;
                    if (j != i)
                        V[j+nt*(l+n*(i+nt*k))] = tmp2;
                    
                }
            }
            
            // End
        }
    }
    Free(lambda);
    Free(S);
    Free(S1);
    Free(U);
    Free(W);
    Free(tmp1);
    Free(elti);
    Free(eltj);
    
}


// REAL
static void simmap_covar_matrix (int *nchar, 
			       double *bt, 
			       double *lambda, 
			       double *S,
                   double *S1,
			       double *sigmasq, 
			       int *nterm, 
			       double *V) {
  double *U, *W;
  double sij, ti, tj;
  double *elti, *eltj;
  double tmp;
  int n = *nchar, nt = *nterm;
  int i, j, k, l, r, s;
  U = Calloc(n*n,double);
  W = Calloc(n*n,double);
  elti = Calloc(n,double);
  eltj = Calloc(n,double);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      U[i+j*n] = 0;
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  U[i+j*n] += S1[i+k*n]*sigmasq[k+l*n]*S1[j+l*n]; // changed S[k+i*n] to S[i+k*n]; to check
	}
      }
    }
  }
  for (i = 0; i < nt; i++) {
    for (j = 0; j <= i; j++) {
      ti = bt[i+i*nt];
      sij = bt[i+j*nt];
      tj = bt[j+j*nt];
      for (k = 0; k < n; k++) {
	elti[k] = exp(-lambda[k]*(ti-sij));
	eltj[k] = exp(-lambda[k]*(tj-sij));
      }
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  V[i+nt*(k+n*(j+nt*l))] = 0; 
	  V[j+nt*(k+n*(i+nt*l))] = 0; 
	  W[k+l*n] = elti[k]*U[k+l*n]*eltj[l]/(lambda[k]+lambda[l]);
        // similar eigenvalues
        // te^t*alpha
	}
      }
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	      tmp = S[k+r*n]*W[r+s*n]*S[l+s*n];
	      V[i+nt*(k+n*(j+nt*l))] += tmp;
	      if (j != i) 
		V[j+nt*(l+n*(i+nt*k))] += tmp;
	    }
	  }
	}
      }
    }
  }
  Free(U);
  Free(W);
  Free(elti);
  Free(eltj);
}

SEXP simmap_covar (SEXP nterm, SEXP bt, SEXP lambda, SEXP S, SEXP S1, SEXP sigmasq) {
    int nprotect = 0;
    SEXP V;
    int nchar, nt, vdim[2];
    nt = INTEGER(nterm)[0];
    nchar = GET_LENGTH(lambda);
    vdim[0] = nt*nchar; vdim[1] = vdim[0];
    PROTECT(V = makearray(2,vdim)); nprotect++;
    
  // Check if there is complex eigenvalues
  if(!isComplex(lambda)){
      
      // REAL
      simmap_covar_matrix(&nchar, REAL(bt), REAL(lambda), REAL(S), REAL(S1), REAL(sigmasq), &nt, REAL(V));
  }else{
     
      // COMPLEX
      simmap_covar_matrix_complex(&nchar, REAL(bt), COMPLEX(lambda), COMPLEX(S), COMPLEX(S1), REAL(sigmasq), &nt, REAL(V));
  }
    
    
  UNPROTECT(nprotect);
  return V;
}


