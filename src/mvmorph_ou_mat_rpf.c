/*-----------Matrice de covariance pour un processus Ornstein-Uhlenbeck-------------------------*/
/*--Matrice stockée au format RPF "column major order" (Fortran Lapack)-------------------------*/
/*-- moins de boucles et calcul plus rapide?----------------------------------------------------*/
/*-mvMORPH 1.0.3 - 2014 - Julien Clavel - julien.clavel@hotmail.fr/julien.clavel@univ-lyon1.fr--*/
#include "mvmorph.h"

// Fixed root covariance matrix
static void mvmorph_covar_OU_RPF_fixed(int *na, double *A, double *ARF, double *alpha, double *sigma){
int i, j, ij, i1, i2, i3, l, n1, nx2, nt, mod, np1x2, n;
double T, sij, ti, tj, tjj, temp, var;
// Paramètres
n = *na;
nt = (1 + n)*n/2; //nbr d'elements dans le format "packed"
mod = n%2;

var=sigma[0]/(2.0*alpha[0]);
    
// taille de la matrice A, N est pair
if(mod == 0){ 
	// Params
	n1 = n / 2;
	np1x2 = n + n + 2;
	ij = nt - n - 1;
    i1 = n1;
	
for (j = n - 1; j >= i1; --j) {
	i2 = j;
	tjj = A[j + j * n];
	
	for (i = 0; i <= i2; ++i) {
	sij = A[i + j * n];
	tj = tjj - sij;
	ti = A[i + i * n] - sij;
	T=ti+tj;
	temp = (1-exp(-2.0*alpha[0]*sij))*exp(-1.0*alpha[0]*T);
	ARF[ij] = temp * var;
	++ij;
	}
	
	i2 = n1 - 1;
	i3 = j - n1;
	tjj = A[i3 + i3 * n];
		
	for (l = j - n1; l <= i2; ++l) {
	sij = A[i3 + l * n];
	ti = A[l + l * n] - sij;
	tj = tjj - sij;
	T=ti+tj;
	temp = (1-exp(-2.0*alpha[0]*sij))*exp(-1.0*alpha[0]*T);
	ARF[ij] = temp * var;
	++ij;
	}
ij -= np1x2;
}// End for j
	
	
// taille de la matrice A, N impair 	
}else{
	// Parameters
	nx2 = n + n;
	n1 = n / 2; // division par un entier (dimension du triangle)
	ij = nt - n;
	i1 = n1;


for (j = n - 1; j >= i1; --j) {
	i2 = j;
	tjj = A[j + j * n];

	for (i = 0; i <= i2; ++i) {
	sij = A[i + j * n];
	tj = tjj - sij;
	ti = A[i + i * n] - sij;
	T=ti+tj;
	temp = (1-exp(-2.0*alpha[0]*sij))*exp(-1.0*alpha[0]*T);
	ARF[ij] = temp * var;
	++ij;
	}
		i2 = n1 - 1;
		i3 = j - n1;
		tjj = A[i3 + i3 * n];
		
	for (l = j - n1; l <= i2; ++l) {
	sij = A[i3 + l * n];
	ti = A[l + l * n] - sij;
	tj =  tjj - sij;
	T=ti+tj;
	temp = (1-exp(-2.0*alpha[0]*sij))*exp(-1.0*alpha[0]*T);
	ARF[ij] = temp * var;
	++ij;
	}
ij -= nx2;
		}
	}// End else
}// End void

// Random root covariance matrix
static void mvmorph_covar_OU_RPF_random(int *na, double *A, double *ARF, double *alpha, double *sigma){
    int i, j, ij, i1, i2, i3, l, n1, nx2, nt, mod, np1x2, n;
    double T, sij, ti, tj, tjj, var;
    // Paramètres
    n = *na;
    nt = (1 + n)*n/2; //nbr d'elements dans le format "packed"
    mod = n%2;
    
    var=sigma[0]/(2.0*alpha[0]);
    
    // taille de la matrice A, N est pair
    if(mod == 0){
        // Params
        n1 = n / 2;
        np1x2 = n + n + 2;
        ij = nt - n - 1;
        i1 = n1;
        
        for (j = n - 1; j >= i1; --j) {
            i2 = j;
            tjj = A[j + j * n];
            
            for (i = 0; i <= i2; ++i) {
                sij = A[i + j * n];
                tj = tjj - sij;
                ti = A[i + i * n] - sij;
                T=ti+tj;
                ARF[ij] = exp(-1.0*alpha[0]*T) * var;
                ++ij;
            }
            
            i2 = n1 - 1;
            i3 = j - n1;
            tjj = A[i3 + i3 * n];
            
            for (l = j - n1; l <= i2; ++l) {
                sij = A[i3 + l * n];
                ti = A[l + l * n] - sij;
                tj = tjj - sij;
                T=ti+tj;
                ARF[ij] = exp(-1.0*alpha[0]*T) * var;
                ++ij;
            }
            ij -= np1x2;
        }// End for j
        
        
        // taille de la matrice A, N impair
    }else{
        // Parameters
        nx2 = n + n;
        n1 = n / 2; // division par un entier (dimension du triangle)
        ij = nt - n;
        i1 = n1;
        
        
        for (j = n - 1; j >= i1; --j) {
            i2 = j;
            tjj = A[j + j * n];
            
            for (i = 0; i <= i2; ++i) {
                sij = A[i + j * n];
                tj = tjj - sij;
                ti = A[i + i * n] - sij;
                T=ti+tj;
                ARF[ij] = exp(-1.0*alpha[0]*T) * var;
                ++ij;
            }
            i2 = n1 - 1;
            i3 = j - n1;
            tjj = A[i3 + i3 * n];
            
            for (l = j - n1; l <= i2; ++l) {
                sij = A[i3 + l * n];
                ti = A[l + l * n] - sij;
                tj =  tjj - sij;
                T=ti+tj;
                ARF[ij] = exp(-1.0*alpha[0]*T) * var;
                ++ij;
            }
            ij -= nx2;
        }
    }// End else
}// End void

// Use random or fixed root covariance matrix in RPF format
SEXP mvmorph_covar_ou_rpf_fixed(SEXP A, SEXP alpha, SEXP sigma) {
 int na; 
	PROTECT(coerceVector(A,REALSXP)); 
	na=INTEGER(GET_DIM(A))[0];
	SEXP ARF;
    PROTECT(ARF = allocVector(REALSXP,(na+1)*na/2)); 
	mvmorph_covar_OU_RPF_fixed(&na,REAL(A),REAL(ARF),REAL(alpha), REAL(sigma));
  UNPROTECT(2);
  return ARF;
}

SEXP mvmorph_covar_ou_rpf_random(SEXP A, SEXP alpha, SEXP sigma) {
    int na;
    PROTECT(coerceVector(A,REALSXP));
    na=INTEGER(GET_DIM(A))[0];
    SEXP ARF;
    PROTECT(ARF = allocVector(REALSXP,(na+1)*na/2));
    mvmorph_covar_OU_RPF_random(&na,REAL(A),REAL(ARF),REAL(alpha), REAL(sigma));
    UNPROTECT(2);
    return ARF;
}

