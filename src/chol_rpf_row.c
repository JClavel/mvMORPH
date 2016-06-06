/* Solve linear system with RPF Cholesky - Julien Clavel - mvMORPH 1.0.3 - 2014 */
/* Fast computation of the Cholesky factorization using BLAS 3 and half memmory */
/* storage (Rectangular Packed Format). The Cholesky factor is the used to      */
/* solve the generalized least square system. "row-major order version"         */

#include "mvmorph.h"


// function to compute the determinant of the matrix A
static void determinant(double *det, double *ARF, int *n){
	int i, j, mod, n1, n2, n3, n4, na=*n;
	#define Br(r) (log(ARF[(r)*n2 + n3 + (r)]))
    #define Br2(r) (log(ARF[(r)*n2 + n4 + (r)]))
    #define Ar(r) (log(ARF[(r)*n1 + (r) + n3]))
    #define Ar2(r) (log(ARF[(r)*n1 + (r) + n4]))	
	// Calcul du determinant
	det[0]=0;
	// Evalue si la matrice est d'ordre pair ou impair
	mod = na%2;
	// taille de la matrice A, N est pair
	if(mod == 0){ 
	n1 = na/2;
	n3 = n1*n1;
	n4 = n3+n1;

		for (i = 0; i < na; i++) {
			if(i < n1){
			det[0]+=Ar(i);
			}else{
			j = i-n1;
			det[0]+=Ar2(j);
			}
		}
        
	det[0]=det[0]*2.0;
	// taille de la matrice A est impair
	}else{ 
	n1 = na/2;
	n2 = na-n1;
	n3 = n1*n2;
	n4 = n3+n2;
		for (i = 0; i < n2; i++) {
			det[0]+=Br(i);
		}
		for(j=0; j < n1; j++){
			det[0]+=Br2(j);
		}
		det[0]=det[0]*2.0;
	}
}

// Add measurement error to the matrix diagonal
static void ms_error(double *A, double *mserr, int *n){
	int nnvcv = *n, i=0, round = down(nnvcv,4);
	#define AA(x) (A[(x)*nnvcv+(x)])
    // loop unrolling
	for(; i<round; i+=4){
		AA(i)+= mserr[i];
		AA(i+1)+= mserr[i+1];
		AA(i+2)+= mserr[i+2];
		AA(i+3)+= mserr[i+3];
	}
	for(; i<nnvcv; i++){
	   AA(i)+= mserr[i];
	   }
}


// Factorisation de Cholesky RPF
SEXP   Chol_RPF(SEXP A, SEXP D, SEXP dat, SEXP nterm, SEXP ndimA, SEXP mserr, SEXP ismserr){
    int n, nt, err, info, one = 1;
    double alpha = 1.;
	char up = 'U', trans = 'T', diag = 'N', side = 'L'; 
	nt = INTEGER(nterm)[0];
	n = INTEGER(ndimA)[0];
	err = INTEGER(ismserr)[0];
	PROTECT(A = coerceVector(A,REALSXP));
	PROTECT(mserr = coerceVector(mserr,REALSXP));
	SEXP Ddat = PROTECT(isReal(dat) ? duplicate(dat): coerceVector(dat, REALSXP));
	SEXP DD = PROTECT(isReal(D) ? duplicate(D): coerceVector(D, REALSXP));
	SEXP ARF = PROTECT(allocVector(REALSXP,(n+1)*n/2));
	SEXP det = PROTECT(allocVector(REALSXP,1));
	
	// add measurement error
	if(err==1){
	ms_error(REAL(A),REAL(mserr), &n);
	}  
	// preparation au format RPF
	F77_CALL(dtrttf)(&trans,&up,&n,REAL(A),&n,REAL(ARF),&info);
    if (info != 0){
        error("the %d argument had an illegal value",info);
    }
    
	// decomposition de Cholesky
	F77_CALL(dpftrf)(&trans,&up,&n, REAL(ARF),&info);
    if (info != 0) {
        if (info > 0) error("the leading minor of order %d is not positive definite",info);
        error("argument %d of Lapack routine %s had invalid value",-info, "dpftrf");
    }

	// systeme lineaire U'x=D
	F77_CALL(dtfsm)(&trans, &side, &up, &trans, &diag, &n, &nt, &alpha, REAL(ARF), REAL(DD), &n);
	// systeme lineaire U'x=dat
	F77_CALL(dtfsm)(&trans, &side, &up, &trans, &diag, &n, &one, &alpha, REAL(ARF), REAL(Ddat), &n);

	// Calcul du determinant
	determinant(REAL(det),REAL(ARF),&n);
	
	// Liste: determinant, cholesky, X
	SEXP vec = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(vec, 0, ARF);
    SET_VECTOR_ELT(vec, 1, det);
	SET_VECTOR_ELT(vec, 2, DD);
	SET_VECTOR_ELT(vec, 3, Ddat);
	
    UNPROTECT (7);
    return (vec); 
}


// Factorisation de Cholesky RPF - avoid explicit computation of U'x=D and U'x=dat
SEXP   Chol_RPF_only(SEXP A, SEXP ndimA, SEXP mserr, SEXP ismserr){
    int n, err, info;
    char up = 'U', trans = 'T';
    
    n = INTEGER(ndimA)[0];
    err = INTEGER(ismserr)[0];
    PROTECT(A = coerceVector(A,REALSXP));
    PROTECT(mserr = coerceVector(mserr,REALSXP));
    SEXP ARF = PROTECT(allocVector(REALSXP,(n+1)*n/2));
    SEXP det = PROTECT(allocVector(REALSXP,1));
    
    // add measurement error
    if(err==1){
        ms_error(REAL(A),REAL(mserr), &n);
    }
    // preparation au format RPF
    F77_CALL(dtrttf)(&trans,&up,&n,REAL(A),&n,REAL(ARF),&info);
    if (info != 0){
        error("the %d argument had an illegal value",info);
    }
    
    // decomposition de Cholesky
    F77_CALL(dpftrf)(&trans,&up,&n,REAL(ARF),&info);
    if (info != 0) {
        if (info > 0) error("the leading minor of order %d is not positive definite",info);
        error("argument %d of Lapack routine %s had invalid value",-info, "dpftrf");
    }
    
    // Calcul du determinant
    determinant(REAL(det),REAL(ARF),&n);
    
    // Liste: cholesky, determinant
    SEXP vec = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(vec, 0, ARF);
    SET_VECTOR_ELT(vec, 1, det);
    
    UNPROTECT (5);
    return (vec);
}

