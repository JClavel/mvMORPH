/*-Matrice de covariance pour un processus Ornstein-Uhlenbeck multivarie--------------*/
/*-mvMORPH 1.0.2 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
#include "covar.h"


static void copy_sparse(int *IA, int *JA, int *nrow, double *V, double *A){

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

// Complex version

static void mvmorph_covar_mat_nult_complex (int *nchar, int *nt, double *bt, Rcomplex *lambda_val, Rcomplex *S_val, double *sigmasq, double *V, Rcomplex *S1_val) {
    double sij, ti, tj, tmp;
    int n = *nchar, nn = *nt;
    int i, j, k, l, s, r, kln, krn, rsn, sln;
    // complex version
    double complex *exp1, *exp2, *exp1l, *exp2l, *G, *U, *F, *lambda, *S, *S1, constzero = 0.0 + 0.0*I;
    U = Calloc(n*n,double complex);
    G = Calloc(n*n,double complex);
    F = Calloc(n*n,double complex);
    exp1l = Calloc(n*n,double complex);
    exp2l = Calloc(n*n,double complex);
    exp1 = Calloc(n*n,double complex);
    exp2 = Calloc(n*n,double complex);
    S = Calloc(n*n,double complex);
    S1 = Calloc(n*n,double complex);
    lambda = Calloc(n,double complex);
    
    for(i=0; i<n; i++){
        lambda[i]=comp(lambda_val[i]);
    }
    
    //zeroing vectors & transform to C complex structure
    for(i = 0; i<n; i++){
        for(j =0; j<n; j++){
            exp1l[i+j*n]= constzero;
            exp2l[i+j*n]= constzero;
            S[i+j*n]=comp(S_val[i+j*n]);
            S1[i+j*n]=comp(S1_val[i+j*n]);
        }
    }
    
    // zeroing the V matrix
    memset(V, 0, ((n*nn)*(n*nn))*sizeof(double));
    
    /* Calcul de la fonction de correlation BBt*/
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            U[i+j*n] = 0;
            for (k = 0; k < n; k++) {
                for (l = 0; l < n; l++) {
                    U[i+j*n] += S1[i+k*n]*sigmasq[k+l*n]*S1[j+l*n]; //S[l+j*n]=>S1[j+l*n]
                }
            }
        }
    }
    
    /* Calcul de la matrice de covariance */
    for (i = 0; i < nn; i++) {
        for (j = 0; j <= i; j++) {
            /* si l'arbre n'est pas ultrametrique */
            sij = bt[j+i*nn];
            ti = bt[i+i*nn]-sij;
            tj = bt[j+j*nn]-bt[i+j*nn];
            
            // Save time: test before
            if(sij!=0.){
                /* preparation du calcul des matrices exponentielles */
                for (k = 0; k < n; k++) {
                    exp1l[k+k*n] = cexp(-lambda[k]*ti);
                    exp2l[k+k*n] = cexp(-lambda[k]*tj);
                }
                
                /* Calcul de G - Gardiner 2004 4.4.47 - Bartoszek et al. 2012 - Meucci  2010 */
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        G[k+l*n] = ((1.0-cexp(-1.0*(lambda[k]+lambda[l])*sij))/(lambda[k]+lambda[l]))*U[k+l*n];
                    }
                }
                
                /* Calcul de SGS - Gardiner 2004 - calcul des matrices exponentielles*/
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        kln = k+l*n;
                        F[kln] = constzero;
                        exp1[kln] = constzero;
                        exp2[kln] = constzero;
                        for (r = 0; r < n; r++) {
                            for (s = 0; s < n; s++) {
                                rsn = r+s*n; krn = k+r*n; sln = s+l*n;
                                F[kln] += S[krn]*G[rsn]*S[l+s*n]; // the Integral
                                exp1[kln] += S[krn]*exp1l[rsn]*S1[sln];
                                exp2[kln] += S[krn]*exp2l[rsn]*S1[sln];
                            }
                        }
                    }
                }
                /* Integrale + matrices exponentielles */
                // Here we can compute directly the product of the real parts as we expect that the matrix exponential of real matrix with complex eigenvalues is real.
                if (j != i) {
                    for (k = 0; k < n; k++) {
                        for (l = 0; l < n; l++) {
                            for (r = 0; r < n; r++) {
                                for (s = 0; s < n; s++) {
                                    tmp = creal(exp1[k+r*n])*creal(F[r+s*n])*creal(exp2[l+s*n]);
                                    V[i+nn*(k+n*(j+nn*l))] += tmp;
                                    V[j+nn*(l+n*(i+nn*k))] += tmp;
                                }
                            }
                        }
                    }
                }else{
                    for (k = 0; k < n; k++) {
                        for (l = 0; l < n; l++) {
                            for (r = 0; r < n; r++) {
                                for (s = 0; s < n; s++) {
                                    tmp = creal(exp1[k+r*n])*creal(F[r+s*n])*creal(exp2[l+s*n]);
                                    V[i+nn*(k+n*(j+nn*l))] += tmp;
                                }
                            }
                        }
                    }
                    //end if diagonal
                }
                //end if zero
            }
            /* Fin de la boucle pour la matrice de covariance
             on libère la mémoire */
        }
    }
    Free(lambda);
    Free(S);
    Free(S1);
    Free(U);
    Free(G);
    Free(F);
    Free(exp1);
    Free(exp2);
    Free(exp1l);
    Free(exp2l);
}


// Real eigenvalues
static void mvmorph_covar_mat_nult (
				   int *nchar, 
                   int *nt,
				   double *bt, 
			       double *lambda, 
			       double *S, 
			       double *sigmasq, 
			       double *V,
                   double *S1) {
   double *U, *G, *F, *exp1, *exp1l, *exp2, *exp2l;
  double sij, ti, tj, tmp;
  int n = *nchar, nn = *nt;
  int i, j, k, l, s, r, n2, n3, n4, kln, krn, rsn, lsn, sn, sln;
  n2 = n*n;
  n3 = n*nn;
  n4 = n3*n3;

  U = Calloc(n2,double);
  G = Calloc(n2,double);
  F = Calloc(n2,double);
  exp1 = Calloc(n2,double);
  exp2 = Calloc(n2,double);
  exp1l = Calloc(n2,double);
  exp2l = Calloc(n2,double);
  //zeroing vectors
memset(exp1l,0,n2);
memset(exp2l,0,n2);
memset(V, 0, n4*sizeof(double));




/* Calcul de la fonction de correlation BBt*/
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      U[i+j*n] = 0;
      for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  U[i+j*n] += S1[i+k*n]*sigmasq[k+l*n]*S1[j+l*n];
	}
      }
    }
  }
  
 // Transformation en Yale matrix sparse format
 // fournir un vecteur de valeurs dans la diagonale
 // fournir le vecteurs d'entrées de la matrice
 // On accede directement aux valeurs dans la matrice et l'on enregistre celles-ci dans le vecteur de résultat @entries slot de l'object as.spam

/* Calcul de la matrice de covariance */
  for (i = 0; i < nn; i++) {
    for (j = 0; j <= i; j++) {
	/* si l'arbre n'est pas ultrametrique */
      sij = bt[j+i*nn];
      ti = bt[i+i*nn]-sij;
      tj = bt[j+j*nn]-sij;
// if no shared history between species we skip calculations to save time
// maybe this if statment is not really efficient... But it save a lot of computations...

	  if(sij!=0.){
	  
/* preparation du calcul des matrices exponentielles */ 
  for (k = 0; k < n; k++) {
	exp1l[k+k*n] = exp(-lambda[k]*ti);
	exp2l[k+k*n] = exp(-lambda[k]*tj);
      }
/* Calcul de G - Gardiner 2004 4.4.47 - Bartoszek et al. 2012 - Najfeld & Havel 1995 */    
 for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  G[k+l*n] = ((1.0-exp(-1.0*(lambda[k]+lambda[l])*sij))/(lambda[k]+lambda[l]))*U[k+l*n];
	}
  }
/* Calcul de SGS - Gardiner 2004 - calcul des matrices exponentielles*/
     for (k = 0; k < n; k++) {
		for (l = 0; l < n; l++) {
		kln = k+l*n;
		       F[kln] = 0;
			exp1[kln] = 0;
		    exp2[kln] = 0;

	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
		krn = k+r*n;
		sn = s*n;
		rsn = r+sn;
		lsn = l+sn;
        sln = s+l*n;
	         F[kln] += S[krn]*G[rsn]*S[lsn];
		  exp1[kln] += S[krn]*exp1l[rsn]*S1[sln];
		  exp2[kln] += S[krn]*exp2l[rsn]*S1[sln];
	    }
	  }
	}
  }
/* Integrale + matrices exponentielles */  

 if (j != i) {
  for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	    tmp = exp1[k+r*n]*F[r+s*n]*exp2[l+s*n];
	    V[i+nn*(k+n*(j+nn*l))] += tmp;
		V[j+nn*(l+n*(i+nn*k))] += tmp;
	    }
	  }
	}
  }
  }else{
    for (k = 0; k < n; k++) {
	for (l = 0; l < n; l++) {
	  for (r = 0; r < n; r++) {
	    for (s = 0; s < n; s++) {
	    tmp = exp1[k+r*n]*F[r+s*n]*exp2[l+s*n];
	    V[i+nn*(k+n*(j+nn*l))] += tmp;
	    }
	  }
	}
  }
  //end if diagonal
  }
  //end if zero
  }
/* Fin de la boucle pour la matrice de covariance
on libère la mémoire */	  
}
}	  

  Free(U);
  Free(G);
  Free(F);
  Free(exp1);
  Free(exp2);
  Free(exp1l);
  Free(exp2l);
}

SEXP mvmorph_covar_ou_sparse (SEXP A, SEXP JA, SEXP IA, SEXP nterm, SEXP bt,SEXP lambda, SEXP S, SEXP sigmasq, SEXP S1) {
 int nprotect = 0;
  SEXP V;
  int nchar, nt, vdim[2], nrow;
  nt = INTEGER(nterm)[0];
  nchar = GET_LENGTH(lambda);
  vdim[0] = nt*nchar; vdim[1] = vdim[0];
  nrow = vdim[0];
  PROTECT(A = coerceVector(A,REALSXP)); nprotect++;
  PROTECT(IA = coerceVector(IA,INTSXP)); nprotect++;
  PROTECT(JA = coerceVector(JA,INTSXP)); nprotect++;
  PROTECT(V = makearray(2,vdim)); nprotect++;
    
    if(!isComplex(lambda)){
        mvmorph_covar_mat_nult(&nchar,&nt,REAL(bt),REAL(lambda),REAL(S),REAL(sigmasq),REAL(V),REAL(S1));
    }else{
        mvmorph_covar_mat_nult_complex(&nchar,&nt,REAL(bt),COMPLEX(lambda),COMPLEX(S),REAL(sigmasq),REAL(V),COMPLEX(S1));
    }

  copy_sparse(INTEGER(IA),INTEGER(JA),&nrow,REAL(V),REAL(A));
  UNPROTECT(nprotect);
  return A;
}
