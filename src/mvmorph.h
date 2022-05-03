/* mvmorph.h 2015-01-01 */
/* Julien Clavel        */
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#define down(x,y) ((x) & ~((y)-1))
#define square(x) (data[(x)]*data[(x)])

extern void F77_CALL(dtrttf)(char *TRANSR, char *UPLO, int *N, double *A,int *LDA, double *ARF, int *INFO);

extern void F77_CALL(dpftrf)(char *TRANSR, char *UPLO, int *N, double *A,int *INFO);

extern void F77_CALL(dtfsm)(char *TRANSR, char *SIDE, char *UPLO, char *TRANS, char *DIAG, int *M, int *N, double *ALPHA, double *A, double *B, int *LDB);

