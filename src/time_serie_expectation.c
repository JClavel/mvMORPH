/*-Expectation for a multivariate Ornstein-Uhlenbeck process time serie---------------*/
/*-mvMORPH 1.0.7 - 2015 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/

#include "functions_complex.h"
#include "covar.h"


// function to compute the matrix exponential assuming a constant expectation for the optimum
static void multi_exp_matrix (int *nvar, int *npoints, double *time, double *lambda, double *S, double *S1, double *matexp) {
    double *expl;
    int n = *nvar, np = *npoints, s_a=*nvar* *nvar;
    int i, j, k, l, m;
    expl = calloc(np*n,sizeof(double));
    
    
    for(m = 0; m < np; m++){
        
        for(i = 0; i<n; i++){
            expl[i+i*n] = exp(-lambda[i]*time[m]);
        }
      
        // Compute the matrix exponential
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                matexp[(k+l*n) + m*s_a] = 0.0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        matexp[(k+l*n) + m*s_a] += S[k+i*n]*expl[i+j*n]*S1[j+l*n];
                    }
                }
            }
        }
        // End
        
    }// end of m
    
    free(expl);
}

// function to compute the optimum expectation through the time serie
static void optimum(int *nvar, int *npoints,const double *time,const double *theta0,const double *theta1, double *expectation,const double *matexp, const double *matdiag){
    
    int i, j, f, n=*nvar, nn=*npoints, s_a=*nvar* *nvar;
    double *res_theta0, *res_theta1, matexp_val;
    
    res_theta0 = calloc(n,sizeof(double));
    res_theta1 = calloc(n,sizeof(double));
    
    // Multiply the vector of optimums by the matrix
    // Multiply the vector of ancestral states by the matrix
    // Sum the exponential
    // i is for the variables (time series)
    // f is for the sampled points
    // j for traversing the matrix A
    
    for(f=0; f<nn; f++){

        for(i=0; i < n; i++){
            //zeroing the temporary vector
            res_theta1[i] = 0.0;
            res_theta0[i] = 0.0;
            
            for(j=0; j < n; j++){
          
            // Matrix-vector product
                matexp_val = matexp[(i+j*n)+f*s_a]; // indice order accession because "rowmajor"
                res_theta0[i] += matexp_val * theta0[j];
                res_theta1[i] += (matdiag[i+j*n]-matexp_val) * theta1[j];
            }
            
             expectation[f+i*nn] = res_theta1[i] + res_theta0[i];
        }
   
    }
    
    free(res_theta0);
    free(res_theta1);
    
}


static void multi_exp_matrix_complex (int *nvar, int *npoints, double *time, Rcomplex *lambda, Rcomplex *S, Rcomplex *S1, double complex *matexp) {
    double complex *expl, *tmp, *tmp2;
    int n = *nvar, np = *npoints, s_a=*nvar* *nvar;
    int i, j, k, l, m, ind, ind1, ind2, ind3, zero=0;
    expl = calloc(np*n,sizeof(double complex));
    tmp = calloc(np*n,sizeof(double complex));
    tmp2 = calloc(np*n,sizeof(double complex));
    
    for(m = 0; m < np; m++){
        
        for(i = 0; i<n; i++){
            
            // complex exponential
            ind=i+i*n;
            cexpti(lambda, expl, time[m], &ind, &i);
            //expl[i+i*n] = cexp(-comp(lambda[i])*time[m]);
            
        }
        
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                ind=(k+l*n) + m*s_a;
                matexp[ind] = 0 + 0*I;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        // 1) complex product
                        ind1=k+i*n; ind2=i+j*n; ind3=j+l*n;
                        cMulti(S, expl, S1, tmp, tmp2, &ind1, &ind2, &ind3, &zero);
                        // 2) sum
                        matexp[ind] += tmp[0];
                    }
                }
            }
        }
        // End
        
    }// end of m
    
    free(expl);
    free(tmp);
    free(tmp2);
}

// function to compute the optimum expectation through the time serie
static void optimum_complex(int *nvar, int *npoints, const double *time, const double *theta0, const double *theta1, double *expectation, double complex *matexp, const double *matdiag){
    
    int i, j, f, n=*nvar, nn=*npoints, s_a=*nvar* *nvar;
    double *res_theta0, *res_theta1, matexp_val;
    
    res_theta0 = calloc(n,sizeof(double));
    res_theta1 = calloc(n,sizeof(double));
    
    // Multiply the vector of optimums by the matrix
    // Multiply the vector of ancestral states by the matrix
    // Sum the exponential
    // i is for the variables (time series)
    // f is for the sampled points
    // j for traversing the matrix A
    
    for(f=0; f<nn; f++){
        //zeroing the temporary vector
        for(i=0; i < n; i++){
            res_theta1[i] = 0.0;
            res_theta0[i] = 0.0;
            
            for(j=0; j < n; j++){
                // extract the real part
                matexp_val = creal(matexp[(i+j*n)+f*s_a]); // take care about the order of eA in "row major" before returning to R
                res_theta0[i] += matexp_val * theta0[j];
                res_theta1[i] += (matdiag[i+j*n]-matexp_val) * theta1[j];
            }
            
            expectation[f+i*nn] = res_theta1[i] + res_theta0[i];
        }
    }
    
    free(res_theta0);
    free(res_theta1);
    
}

SEXP Expect_matrix(SEXP S1, SEXP S, SEXP lambda, SEXP time, SEXP theta0, SEXP theta1, SEXP matdiag){
    
    int nvar, npoints, nprotect=0;
    
    nvar = GET_LENGTH(lambda);
    npoints = GET_LENGTH(time);
    
    PROTECT(time = coerceVector(time,REALSXP)); nprotect++;
    PROTECT(theta0 = coerceVector(theta0,REALSXP)); nprotect++;
    PROTECT(theta1 = coerceVector(theta1,REALSXP)); nprotect++;
    // results
    SEXP expectation = PROTECT(allocVector(REALSXP,nvar*npoints)); nprotect++;
    
    if(!isComplex(lambda)){
    // eigenvectors
    PROTECT(S1 = coerceVector(S1,REALSXP)); nprotect++;
    PROTECT(S = coerceVector(S,REALSXP)); nprotect++;
    // matrix exponential
    SEXP matexp = PROTECT(allocVector(REALSXP,nvar*nvar*npoints)); nprotect++;

    // Compute the exponential matrix
    multi_exp_matrix (&nvar, &npoints, REAL(time), REAL(lambda), REAL(S), REAL(S1), REAL(matexp));
    
    // Compute the expectations
    optimum (&nvar, &npoints, REAL(time), REAL(theta0), REAL(theta1), REAL(expectation), REAL(matexp), REAL(matdiag));
    // Done.
        
    }else{
    
    double complex *matexp;
    // complex eigenvalues & eigenvectors
    PROTECT(S1 = coerceVector(S1,CPLXSXP)); nprotect++;
    PROTECT(S = coerceVector(S,CPLXSXP)); nprotect++;
        
    // alloc a complex vector in C rather than R structure...
    matexp = calloc(nvar*nvar*npoints,sizeof(double complex));
        
    // Compute the exponential matrix
    multi_exp_matrix_complex (&nvar, &npoints, REAL(time), COMPLEX(lambda), COMPLEX(S), COMPLEX(S1), matexp);
        
    // Compute the expectations
    optimum_complex(&nvar, &npoints, REAL(time), REAL(theta0), REAL(theta1), REAL(expectation), matexp, REAL(matdiag));
    // Done.
    // Free the memory
    free(matexp);
    }


    UNPROTECT(nprotect);
    return expectation;
    
}


// Matrix of Weight for the multivariate OU (REAL eigenvalues)
static void build_w(int *nvar, int *npoints, const double *time, double *expectation, const double *matexp, const double *matdiag){
    
    int i, j, f, k, n=*nvar, nn=*npoints, s_a=*nvar* *nvar, s_b=*nvar* *npoints*2, s_c=*nvar* *npoints;
    double matexp_val, matexp_val1;
    
    for(f=0; f<nn; f++){
        
        for(i=0, k=0 ; i < n; i++, k+=nn){
            
            for(j=0; j < n; j++){
                
                // Matrix-vector product
                matexp_val = matexp[(i+j*n)+f*s_a]; // indice order accession because "rowmajor"
                matexp_val1 = (matdiag[i+j*n]-matexp_val);
                
                // Make the W matrix, for theta 0 and 1.
                expectation[k+j*s_b+f] = matexp_val;
                expectation[k+j*s_b+f+s_c] = matexp_val1;
            }
    
        }
        
    }
 
}


// Matrix of Weight for the multivariate OU
static void build_w_complex(int *nvar, int *npoints, const double *time, double *expectation, const double complex *matexp, const double *matdiag){
    
    int i, j, f, k, n=*nvar, nn=*npoints, s_a=*nvar* *nvar, s_b=*nvar* *npoints*2, s_c=*nvar* *npoints;
    double matexp_val, matexp_val1;
    
    for(f=0; f<nn; f++){
        
        for(i=0, k=0 ; i < n; i++, k+=nn){
            
            for(j=0; j < n; j++){
                
                // Matrix-vector product
                matexp_val = creal(matexp[(i+j*n)+f*s_a]); // indice order accession because "rowmajor"
                matexp_val1 = (matdiag[i+j*n]-matexp_val);
                
                // Make the W matrix, for theta 0 and 1.
                expectation[k+j*s_b+f] = matexp_val;
                expectation[k+j*s_b+f+s_c] = matexp_val1;
            }
            
        }
        
    }
    
}

// Weight matrix
SEXP Weight_matrix(SEXP S1, SEXP S, SEXP lambda, SEXP time, SEXP matdiag){
    
    int nvar, npoints, vdim[2], nprotect = 0;
    nvar = GET_LENGTH(lambda);
    npoints = GET_LENGTH(time);

    SEXP expectation;
    vdim[0] = npoints*nvar; vdim[1] = nvar*2;
    PROTECT(expectation = makearray(2,vdim)); nprotect++;
    
    if(!isComplex(lambda)){
        // eigenvectors
        PROTECT(S1 = coerceVector(S1,REALSXP)); nprotect++;
        PROTECT(S = coerceVector(S,REALSXP)); nprotect++;
        // matrix exponential
        SEXP matexp = PROTECT(allocVector(REALSXP,nvar*nvar*npoints)); nprotect++;
    
        // Compute the exponential matrix
        multi_exp_matrix (&nvar, &npoints, REAL(time), REAL(lambda), REAL(S), REAL(S1), REAL(matexp));
    
        // Compute the expectations
        build_w (&nvar, &npoints, REAL(time), REAL(expectation), REAL(matexp), REAL(matdiag));
        // Done.

    }else{
        
        double complex *matexp;
        // complex eigenvalues & eigenvectors
        PROTECT(S1 = coerceVector(S1,CPLXSXP)); nprotect++;
        PROTECT(S = coerceVector(S,CPLXSXP)); nprotect++;
        
        // alloc a complex vector in C rather than R structure...
        matexp = calloc(nvar*nvar*npoints,sizeof(double complex));
        
        // Compute the exponential matrix
        multi_exp_matrix_complex (&nvar, &npoints, REAL(time), COMPLEX(lambda), COMPLEX(S), COMPLEX(S1), matexp);
        
        // Compute the expectations
        build_w_complex (&nvar, &npoints, REAL(time), REAL(expectation), matexp, REAL(matdiag));
        
        // Done.
        // Free the memory
        free(matexp);
    }


    UNPROTECT(nprotect);
    return expectation;
    
}






