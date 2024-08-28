/* Solve linear system with RPF Cholesky - Julien Clavel - mvMORPH 1.0.3 - 2014 */
/* Fast computation of the quadratic product e'V^-1e with the residuals and the */
/* factorized matrix A.                                                         */
/* The result is used to compute the log-likelihood                             */

#include "mvmorph.h"

// function to compute the quadratic product
SEXP   Chol_RPF_quadprod(SEXP U, SEXP resid, SEXP nterm){
    int n, info = 0, one = 1;
    double alpha = 1.;
	char up = 'U', trans = 'T', diag = 'N', side = 'L'; 
	n = INTEGER(nterm)[0];
	PROTECT(U = coerceVector(U,REALSXP));
	SEXP Ddat = PROTECT(isReal(resid) ? duplicate(resid): coerceVector(resid, REALSXP));
	SEXP Bet = PROTECT(allocVector(REALSXP,1));
    double *beta = REAL(Bet), *data = REAL(Ddat), *chol = REAL(U);
   // systeme lineaire U'x=dat
	F77_CALL(dtfsm)(&trans, &side, &up, &trans, &diag, &n, &one, &alpha, chol, data, &n FCONE FCONE FCONE FCONE FCONE);
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
