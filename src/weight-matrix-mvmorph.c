/*-Matrice de poids W pour un processus Ornstein-Uhlenbeck multivarie-----------------*/
/*-mvMORPH 1.0.3 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/
/*-scale the weight matrix so that the rows sum to 1----------------------------------*/
/*-generalization to non-ultrametric trees and non-symmetric A matrix-----------------*/
/*-Parts of codes modified with permission from A. King OUCH -------------------------*/

#include "covar.h"
#include "functions_complex.h"


// Complex case
static void weight_matrix_complex (int *nchar, int *neps, double *epochs, Rcomplex *lambda, Rcomplex *S, Rcomplex *S1, double complex *y) {
    double complex *elt, *tmp;
    double t;
    int n = *nchar, np = *neps;
    int i, j, k, r, ind, ind1, ind2, ind3, zero=0;
    elt = calloc(n*np,sizeof(double complex));
    tmp = calloc(1,sizeof(double complex));
    
    // Prepare the exponentials for lambda
    for (i = 0; i < np; i++) {
        t = epochs[0]-epochs[i];
        for (j = 0; j < n; j++){
            // complex exponential
            ind1=i+np*j;
            cexpti(lambda, elt, t, &ind1, &j);
            // elt[i+np*j] = cexp(-comp(lambda[j])*t);
        }
    }
    
    // substract the exponentials
    for (i = 0; i < np-1; i++) {
        for (j = 0; j < n; j++) {
            elt[i+np*j] -= elt[i+1+np*j];
        }
    }
    
    // compute the matrix exponential PDPt and save it in y
    for (i = 0; i < np; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                ind=j+n*(k+n*i);
                y[ind] = 0.0 + 0.0*I;
                
                for (r = 0; r < n; r++){
                    // 1) mult (S*elt)*S1
                    ind1=j+n*r; ind2=i+np*r; ind3=r+n*k;
                    cMulti(S, elt, S1, tmp, tmp, &ind1, &ind2, &ind3, &zero);
                    // 2) sum to y
                    y[ind] += tmp[0];
                }
            }
        }
    }
    
    free(elt);
    free(tmp);
}



// Real case

static void multi_weight_matrix (int *nchar, int *neps, double *epochs, double *lambda, double *S, double *S1, double *y) {
  double *elt;
  double t;
  int n = *nchar, np = *neps;
  int i, j, k, r;
  elt = calloc(np*n,sizeof(double));
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
    
  free(elt);
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
  double *wp, *bp;
  int nchar, nt, *nreg, totreg, np, xdim[2], ptr, thetaO;
  int i, j, k, n, q;
  nchar = GET_LENGTH(lambda);  
  nt = INTEGER(nterm)[0];
  thetaO = INTEGER(root)[0];
  nreg = calloc(nchar,sizeof(int));
  totreg = 0;
    
    for (i = 0; i < nchar; i++) {
        nreg[i] = INTEGER(GET_DIM(VECTOR_ELT(VECTOR_ELT(beta,0),i)))[1];
        totreg += nreg[i];
    }
    
        xdim[0] = nt*nchar; xdim[1] = totreg;
        PROTECT(W = makearray(2,xdim)); nprotect++;
    
    if(!isComplex(lambda)){
        
        double *y;
        
        for (i = 0; i < nt; i++) {
            np = GET_LENGTH(VECTOR_ELT(epochs,i));
            y = calloc(nchar*nchar*np,sizeof(double));
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
            
            
            free(y);
        }
    }else{
        double complex *y;
        
        for (i = 0; i < nt; i++) {
            np = GET_LENGTH(VECTOR_ELT(epochs,i));
            // alloc a dynamic complex vector...
            y = calloc(nchar*nchar*np,sizeof(double complex));
            weight_matrix_complex(&nchar,&np,REAL(VECTOR_ELT(epochs,i)),COMPLEX(lambda), COMPLEX(S), COMPLEX(S1), y);
            for (n = 0, ptr = 0; n < nchar; ptr += nt*nchar*nreg[n++]) {
                wp = &(REAL(W))[ptr];
                bp = REAL(VECTOR_ELT(VECTOR_ELT(beta,i),n));
                for (j = 0; j < nchar; j++) {
                    for (k = 0; k < nreg[n]; k++) {
                        wp[i+nt*(j+nchar*k)] = 0;
                        for (q = 0; q < np; q++) {
                            // extract only the real part of the matrix exponential
                            wp[i+nt*(j+nchar*k)] += creal(y[n+nchar*(j+nchar*q)])*bp[q+np*k];
                        }
                    }
                }
            }
            free(y);
        }

    }
        free(nreg);
    
    // Row standardize the matrix if the root is not estimated
    if(thetaO==1){
        int nsp;
        nsp=nt*nchar;
        row_stand(REAL(W),&nsp,&totreg,&nsp);
    }
    
  UNPROTECT(nprotect);
  return W;
}

