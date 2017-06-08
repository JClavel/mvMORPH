// Prunning algorithm to compute the matrix square root of a phylogenetic tree covariance
// See details in Stone 2011 - Syst. Biol.; and also Khabbazian et al. 2016 - Meth. Ecol. Evol.
// Julien Clavel - 2017; clavel@biologie.ens.fr
// NOTE : we need the transpose of D to decorrelate the traits
// Return the matrix square root, the variance at the root, and variance at each nodes

#include <R.h>
#include <Rinternals.h>

static void phylo_squareRoot(const int *nsp, const int *edge1, const int *edge2, double *tempbl, double *F, double *D, double *var_prun, double *root_v, double *V, int *invMat){
    int i, j, k, l, f, anc, da, d1, d2, ntip, indice;
    double sumbl, tfinal, t1, t2;
    
    ntip=*nsp;
    indice = 0; // counter for the columns of D (the matrix square root)
    
    // Matrix square root = chol(C)
    if(*invMat==0){
        
        for (i = 0; i < ntip * 2 - 3; i += 2) {
            /*
            ntip*2-3 is the dim of the edge list. Because we take j to be i+1, we must loop up to nedge-1
            D is ntip*ntip
            F is ntip*2-1*ntip; with first part of the structure a diagonal matrix of size ntip*ntip
            */
            j = i + 1;
        
            anc = edge1[i];                     // ancestor
            da = anc - 1;
            d1 = edge2[i] - 1;                  // first descent
            d2 = edge2[j] - 1;                  // 2nd descent
            t1 = tempbl[i];                     // br length for d1
            t2 = tempbl[j];                     // br length for d2
            sumbl = t1 + t2;                    // total br length
            var_prun[anc - ntip - 1] = sumbl;   // variance of "contrasts"
        
            // update the matrix
            for(f=0; f<ntip; f++){
            
                // Square root
                D[indice*ntip + f] = (F[d1*ntip + f]*t1 - F[d2*ntip + f]*t2)/sqrt(sumbl);
            
                // update the matrix F;
                F[da*ntip + f] = F[d1*ntip + f] + F[d2*ntip + f];
            
            }
        
            /* find the edge where `anc' is a descendant (except if at the root):
            it is obviously below the j'th edge */
            if (j != ntip * 2 - 3) {
                k = j + 1;
                while (edge2[k] != anc) k++;
                tempbl[k] += t1*t2/sumbl;
            }
        
            indice++;
        }
    
        // The last column of the square root matrix
        tfinal = 1./(1./t1 + 1./t2);
        for(f=0; f<ntip; f++){
            D[indice*ntip + f] = F[(1+indice)*ntip + f]*sqrt(tfinal);
        }
        
    // Matrix square root of the inverse = chol(C^-1)
    }else{
        
        
        for (i = 0; i < ntip * 2 - 3; i += 2) {
            /*
             ntip*2-3 is the dim of the edge list. Because we take j to be i+1, we must loop up to nedge-1
             D is ntip*ntip
             F is ntip*2-1*ntip; with first part of the structure a diagonal matrix of size ntip*ntip
             */
            j = i + 1;
            
            anc = edge1[i];                     // ancestor
            da = anc - 1;
            d1 = edge2[i] - 1;                  // first descent
            d2 = edge2[j] - 1;                  // 2nd descent
            t1 = tempbl[i];                     // br length for d1
            t2 = tempbl[j];                     // br length for d2
            sumbl = t1 + t2;                    // total br length
            var_prun[anc - ntip - 1] = sumbl;   // variance of "contrasts"
            
            // update the matrix
            for(f=0; f<ntip; f++){
                
                // Square root inverse
                D[indice*ntip + f] = (F[d1*ntip + f] - F[d2*ntip + f])/sqrt(sumbl);
                
                // update the matrix F; we can use it in the same loop as we are updating a column different from d1 and d2?
                F[da*ntip + f] = (F[d1*ntip + f]*t2 + F[d2*ntip + f]*t1)/sumbl;
                
            }
            
            /* find the edge where `anc' is a descendant (except if at the root):
             it is obviously below the j'th edge */
            if (j != ntip * 2 - 3) {
                k = j + 1;
                while (edge2[k] != anc) k++;
                tempbl[k] += t1*t2/sumbl;
            }
            
            indice++;
        }
        
        // The last column of the square root matrix
        tfinal = 1./(1./t1 + 1./t2);
        for(f=0; f<ntip; f++){
            D[indice*ntip + f] = F[(1+indice)*ntip + f]/sqrt(tfinal);
        }
    }
    
    // End
    // NOTE: we compute the variance at the root
    k=0;
    l=0;
    while(k !=2){
        if(edge1[l] == ntip+1){
            root_v[k]=tempbl[l];
            k++;
        }
        l++;
    }
    
    V[0]= root_v[0]*root_v[1]/(root_v[0]+root_v[1]);
    
}


/* Main function to compute the variance terms and matrix square root*/

SEXP squareRootM(SEXP edge1, SEXP edge2, SEXP edgelength, SEXP nsp, SEXP inverse){
    int ntip, i, dim2;

    ntip = INTEGER(nsp)[0];
    dim2 = (ntip*2-1);
    
    // Protect vectors
    PROTECT(edge1 = coerceVector(edge1,INTSXP));
    PROTECT(edge2 = coerceVector(edge2,INTSXP));

    // copy the branch lengths
    SEXP tempblength = PROTECT(isReal(edgelength) ? duplicate(edgelength): coerceVector(edgelength, REALSXP));
    
    // make the matrix to store the results
    SEXP V1 = PROTECT(allocVector(REALSXP,ntip-1));
    SEXP V0 = PROTECT(allocVector(REALSXP,1));
    SEXP root_v = PROTECT(allocVector(REALSXP,2));
    SEXP D = PROTECT(allocMatrix(REALSXP,ntip,ntip));
    SEXP F = PROTECT(allocVector(REALSXP,dim2*ntip));
    memset(REAL(D), 0, (ntip*ntip) * sizeof(double));
    memset(REAL(F), 0, (ntip*dim2) * sizeof(double));
    
    // Fill the diagonal part of F
    double *ident = REAL(F);
    for(i=0; i<ntip; i++) ident[i*ntip + i] = 1;
    
    // Compute the square root
      phylo_squareRoot(&ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(F), REAL(D), REAL(V1), REAL(root_v), REAL(V0), INTEGER(inverse));
    
    // resultats
    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result, 0, D);
    SET_VECTOR_ELT(result, 1, V1);
    SET_VECTOR_ELT(result, 2, V0);
    
    UNPROTECT(9);
    return (result);
    
}
