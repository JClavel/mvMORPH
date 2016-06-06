/* functions_complex.h 2016-01-01 */
/* Julien Clavel - mvMORPH        */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <complex.h>
#define sq(x) ((x)*(x))
// transform complex R structure to C structure
#define comp(x) ((x.r) + ((x.i)*I))

// complex exponential of negative eigenvalues
static void cexpti(Rcomplex *eigvalues, double complex *elt, double t, int *ind, int *ind2){
    int n=*ind, i=*ind2;
    // exponential with complex numbers Gockenbach 2010, p. 259
    elt[n] = exp(-eigvalues[i].r*t) * (cos(-eigvalues[i].i*t) + sin(-eigvalues[i].i*t)*I);
}

// Mix the C and R complex objects... allow the use of dynamic vectors. To optimize later with special structures.
static void cMulti(Rcomplex *mat1, double complex *mat2, Rcomplex *mat3, double complex *results, double complex *tmp, int *ind1, int *ind2, int *ind3, int *ind4){
    tmp[0] = comp(mat1[*ind1]) * mat2[*ind2];
    results[*ind4] = tmp[0] * comp(mat3[*ind3]);
}
