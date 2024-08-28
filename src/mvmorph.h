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

La_extern void F77_NAME(dtrttf)(const char* transr, const char* uplo, const int* n,
         const double* a, const int* lda,
         double* arf, int* info FCLEN FCLEN);

La_extern void F77_NAME(dpftrf)(const char* transr, const char* uplo, const int* n,
         double* a, int* info FCLEN FCLEN);

La_extern void F77_NAME(dtfsm)(const char* transr, const char* side, const char* uplo, const char* trans, const char* diag,
        const int* m, const int* n, const double* alpha, const double* a,
        double* b, const int* ldb FCLEN FCLEN FCLEN FCLEN FCLEN);


