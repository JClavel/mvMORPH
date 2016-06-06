/*-Likelihood computation with the contrast method for multivariate models------------*/
/*-mvMORPH 1.0.3 - 2014 - Julien Clavel - julien.clavel@hotmail.fr--------------------*/

#include "mvmorph.h"

// Copy of a vector 
static void copybrlength(int *ntip, double *blength, double *tempblength){
#define down(x,y) ((x) & ~((y)-1))
int veclength=*ntip*2-2;
int i=0, round=down(veclength,4);
	
    for(; i<round; i+=4){
        tempblength[i]=blength[i];
        tempblength[i+1]=blength[i+1];
        tempblength[i+2]=blength[i+2];
        tempblength[i+3]=blength[i+3];
   }
    for(;i<veclength; i++){
        tempblength[i]=blength[i];
    }
   	
}

// Copy phenotype
static void copypheno(int *ntip, int *ntraits, int *ntotal, const double *pheno, double *phenocopy){
    int i, j, nsp=*ntip, ntr=*ntraits, ntot=*ntotal ;
    
    // Préparation du jeu de données
    for(i=0; i<ntr; i++){
        for(j=0; j<nsp; j++){
            phenocopy[j+i*ntot]=pheno[j+i*nsp];
        }
    }
}

// Early Burst tree transformation
static void ebTree(double *times, double *tempblength, double *rate, int *edge1, int *edge2, int *ntip){
    int i, nt=*ntip, ind, a1, d1;
    ind=nt*2-2;
    
    for(i=0; i<ind; i++){
        d1=edge2[i]-1;
        a1=edge1[i]-1;
        tempblength[i]=(exp(rate[0]*times[d1])-exp(rate[0]*times[a1]))/rate[0];
        
    }
}

// Ornstein-Uhlenbeck tree transformation
static void ouTree(double *times, double *tempblength, double *blength, double *Tmax, double *rate, int *edge1, int *ntip){
    double bl, age, t1, t2, var;
    int inded1, ind, nt=*ntip, i;
    
    ind=nt*2-2;
    var=1./(2.*rate[0]);
    
    for(i=0; i<ind; i++){
        bl=blength[i];
        inded1=edge1[i];
        age=times[inded1-(nt+1)];
        t1=Tmax[0]-age;
        t2=t1+bl;
        tempblength[i]=var*exp(-2.*rate[0] * (Tmax[0]-t2)) * (1. - exp(-2. * rate[0] * t2))-var*exp(-2.*rate[0] * (Tmax[0]-t1)) * (1. - exp(-2. * rate[0] * t1));
    }
}

// quadratic product
static void dotprodX(double *S, double *iZ, double *contr, int *ntraits, int *numbnod){
int i, j, f, n=*ntraits, nnod=*numbnod;
double *res;

res = Calloc(n,double);

for(f=0; f<nnod; f++){
	//zeroing the temporary vector
	for(i=0; i < n; i++){
	res[i] = 0.0;
		for(j=0; j < n; j++){
	res[i] += iZ[j+i*n] * contr[f+j*nnod];
		}
	}
	for(i=0; i<n; i++){
		S[0]+=res[i]*contr[f+i*nnod];
	}
}
Free(res);
}

 // calcul des contrastes pour chacun des traits (modifié d'après "pic" dans APE)
  /* pic.c       2006-11-13 */
  /* Copyright 2006 Emmanuel Paradis */
  /* This file is part of the R-package `ape'. */
  /* See the file ../COPYING for licensing issues. */

static void phylo_pic(int *ind, int *ntotal, int *numbnode, int *nsp, int *edge1, int *edge2, double *tempblength, double *pheno, double *var_contr, double *ancstates, double *root_v, double *V, double *contr){
    int i, j, k, l, ic, anc, d1, d2, ntot, numbnod, ntip, f;
    double sumbl;
    ntot=*ntotal;
    numbnod=*numbnode;
    ntip=*nsp;
    f=*ind;
    
    for (i = 0; i < ntip * 2 - 3; i += 2) {
        j = i + 1;
        anc = edge1[i];
        d1 = edge2[i] - 1;
        d2 = edge2[j] - 1;
        sumbl = tempblength[i] + tempblength[j];
        ic = anc - ntip - 1;
        contr[ic+f*numbnod] = (pheno[d1+f*ntot] - pheno[d2+f*ntot])/sqrt(sumbl);
        var_contr[ic+f*numbnod] = sumbl;
        pheno[(anc - 1)+f*ntot] = (pheno[d1+f*ntot]*tempblength[j] + pheno[d2+f*ntot]*tempblength[i])/sumbl;
        /* find the edge where `anc' is a descendant (except if at the root):
         it is obviously below the j'th edge */
        if (j != ntip * 2 - 3) {
            k = j + 1;
            while (edge2[k] != anc) k++;
            tempblength[k] += tempblength[i]*tempblength[j]/sumbl;
        }
    }
    
    // Calcul de l'etat ancestral
    ancstates[f]=pheno[ntip+f*ntot];
    // Calcul de la Variance
    k=0;
    l=0;
    while(k !=2){
        if(edge1[l] == ntip+1){
            root_v[k]=tempblength[l];
            k++;
        }
        l++;
    }
    
         V[f]= root_v[0]*root_v[1]/(root_v[0]+root_v[1]);
}

static void phylo_pic_single(int *ind, int *ntotal, int *numbnode, int *nsp, int *edge1, int *edge2, double *tempblength, double *pheno, double *var_contr, double *ancstates, double *root_v, double *V, double *contr){
    int i, j, k, l, ic, anc, d1, d2, ntot, numbnod, ntip, f, ntraits;
    double sumbl;
    ntot=*ntotal;
    numbnod=*numbnode;
    ntip=*nsp;
    ntraits=*ind;
    
    for (i = 0; i < ntip * 2 - 3; i += 2) {
        j = i + 1;
        anc = edge1[i];
        d1 = edge2[i] - 1;
        d2 = edge2[j] - 1;
        sumbl = tempblength[i] + tempblength[j];
        ic = anc - ntip - 1;
        
        for(f=0; f<ntraits; f++){
            contr[ic+f*numbnod] = (pheno[d1+f*ntot] - pheno[d2+f*ntot])/sqrt(sumbl);
            var_contr[ic+f*numbnod] = sumbl;
            pheno[(anc - 1)+f*ntot] = (pheno[d1+f*ntot]*tempblength[j] + pheno[d2+f*ntot]*tempblength[i])/sumbl;
        }
        /* find the edge where `anc' is a descendant (except if at the root):
         it is obviously below the j'th edge */
        if (j != ntip * 2 - 3) {
            k = j + 1;
            while (edge2[k] != anc) k++;
            tempblength[k] += tempblength[i]*tempblength[j]/sumbl;
        }
    }
    
    // Calcul de la Variance
    k=0;
    l=0;
    while(k !=2){
        if(edge1[l] == ntip+1){
            root_v[k]=tempblength[l];
            k++;
        }
        l++;
    }
    
    for(f=0 ; f<ntraits; f++){
        V[f]= root_v[0]*root_v[1]/(root_v[0]+root_v[1]);
    
        // Calcul de l'etat ancestral
        ancstates[f]=pheno[ntip+f*ntot];
    }
}

// Models for Phylogenetic Independent Contrasts likelihood
// 1-Early Burst/ACDC model
// 2-Ornstein-Uhlenbeck process
// 3-Brownian Motion
// 4-Generalized Brownian Motion
// 5-Generalized Brownian Motion - ancestral state estimation
// 6-Generalized Brownian Motion on a single topology
// 7-Generalized Brownian Motion on a single topology- ancestral state estimation
// 8-Generalized Brownian Motion on a single topology- sigma estimation
// 9-Generalized Brownian Motion - sigma estimation
// 10-OU process with user mean
// 11-default - Brownian Motion on a single topology


SEXP PIC_gen(SEXP x, SEXP n, SEXP Nnode, SEXP nsp, SEXP edge1, SEXP edge2, SEXP edgelength, SEXP times, SEXP rate, SEXP Tmax, SEXP Model, SEXP mu, SEXP sigma){
  int nodnbtr, numbnod, ntip, i, ntot, ntraits, dimrtrait, f, model, info = 0, neg = 0;
  char transa = 'T', transb = 'N';
  double one = 1.0, zero = 0.0;
  numbnod = INTEGER(Nnode)[0];
  ntip = INTEGER(nsp)[0];
  ntraits = INTEGER(n)[0];
  model = INTEGER(Model)[0];
  ntot = ntip+numbnod;
  nodnbtr = numbnod*ntraits;
  
  //Allocation des vecteurs
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(times = coerceVector(times,REALSXP));
  PROTECT(edge1 = coerceVector(edge1,INTSXP));
  PROTECT(edge2 = coerceVector(edge2,INTSXP));
  PROTECT(edgelength = coerceVector(edgelength,VECSXP));
  SEXP tempblength = PROTECT(allocVector(REALSXP,ntip*2-2));
  SEXP contr = PROTECT(allocVector(REALSXP,nodnbtr));  
  SEXP var_contr = PROTECT(allocVector(REALSXP,nodnbtr));  
  SEXP pheno = PROTECT(allocVector(REALSXP,ntot*ntraits));
  SEXP ancstates = PROTECT(allocVector(REALSXP,ntraits));
  SEXP muRoot = PROTECT(allocVector(REALSXP,ntraits));
  SEXP V = PROTECT(allocVector(REALSXP,ntraits));
  SEXP root_v = PROTECT(allocVector(REALSXP,2));
  SEXP var = PROTECT(allocVector(REALSXP,1));
  memset(REAL(pheno),0,(ntot)*sizeof(double));

 
  // Préparation du jeu de données
  copypheno(&ntip, &ntraits, &ntot, REAL(x), REAL(pheno));
     
   // arbre transformé par traits
   //switch depending on branch transformation
     switch (model) {
         case 1:
             for(f=0; f < ntraits; f++){
                 if(REAL(rate)[f]==0.0) model=11;
             ebTree(REAL(times),REAL(tempblength),REAL(rate),INTEGER(edge1), INTEGER(edge2), &ntip);
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
             }
             break;
             
         case 2:
             for(f=0; f < ntraits; f++){
                 if(REAL(rate)[f]==0.0) model=11;
             ouTree(REAL(times),REAL(tempblength),REAL(VECTOR_ELT(edgelength,0)),REAL(Tmax),REAL(rate),INTEGER(edge1), &ntip);
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
             }
             break;
             
         case 3:
            for(f=0; f < ntraits; f++){
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,f)), REAL(tempblength)); // estimate sigma & theta
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
            }
             break;
             
         case 4:
            for(f=0; f < ntraits; f++){
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,f)), REAL(tempblength)); // user sigma & theta
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
            }
             break;
             
         case 5:
            for(f=0; f < ntraits; f++){
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,f)), REAL(tempblength)); // user sigma estimated theta
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
            }
             break;
             
         case 6:
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,0)), REAL(tempblength)); // 1 topo. user sigma & theta
             phylo_pic_single(&ntraits, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
             break;
             
         case 7:
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,0)), REAL(tempblength)); // 1 topo. user sigma, estimated theta
             phylo_pic_single(&ntraits, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
             break;
             
         case 8:
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,0)), REAL(tempblength)); // 1 topo. estimated sigma, user theta
             phylo_pic_single(&ntraits, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
             break;
             
         case 9:
            for(f=0; f < ntraits; f++){
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,f)), REAL(tempblength)); // estimated sigma, user theta
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
            }
             break;
             
         case 10:
            for(f=0; f < ntraits; f++){
                 if(REAL(rate)[f]==0.0) model=11;
             ouTree(REAL(times),REAL(tempblength),REAL(VECTOR_ELT(edgelength,0)),REAL(Tmax),REAL(rate),INTEGER(edge1), &ntip);
             phylo_pic(&f, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
            }
             break;

             
         default:
             copybrlength(&ntip, REAL(VECTOR_ELT(edgelength,0)), REAL(tempblength));
             phylo_pic_single(&ntraits, &ntot, &numbnod, &ntip, INTEGER(edge1), INTEGER(edge2), REAL(tempblength), REAL(pheno), REAL(var_contr), REAL(ancstates), REAL(root_v), REAL(V), REAL(contr));
     }
 
    
 
    /* Compute the likelihood */
	
	dimrtrait = ntraits*ntraits;
	//SEXP Z = PROTECT(allocVector(REALSXP,dimrtrait));
    SEXP Z = PROTECT(allocMatrix(REALSXP,ntraits,ntraits));
    
    // Sigma matrix is provided
    if(model==1 || model==2 || model==4 || model==5 || model==6 || model==7 || model==10){
        for(i=0; i<dimrtrait; i++) REAL(Z)[i]=REAL(sigma)[i];
    }else{
    // crossproduct of contrast matrix
	F77_CALL(dgemm)(&transa, &transb, &ntraits, &ntraits, &numbnod, &one, REAL(contr), &numbnod, REAL(contr), &numbnod, &zero, REAL(Z), &ntraits);
	// Division par le nombre de taxons (ou taxon -1 pour REML)
		for(i =0; i<dimrtrait; i++) REAL(Z)[i] /= ntip;
    }
	// LU factorization
	SEXP iZ = PROTECT(isReal(Z) ? duplicate(Z): coerceVector(Z, REALSXP));
	SEXP IPIV = PROTECT(allocVector(INTSXP,ntraits));
	F77_CALL(dgetrf)(&ntraits, &ntraits, REAL(iZ), &ntraits, INTEGER(IPIV), &info);
	// Computation of the determinant
	SEXP det = PROTECT(allocVector(REALSXP,1));
	REAL(det)[0]=1.0;
	for(i=0; i<ntraits; i++){
	REAL(det)[0]*=REAL(iZ)[i+i*ntraits];
	 if (INTEGER(IPIV)[i] != (i+1)) neg = !neg;
	}
	// the sign of the determinant
	REAL(det)[0] = neg? -REAL(det)[0]:REAL(det)[0];
	// log determinant
	REAL(det)[0]=log(REAL(det)[0]);
	
	// Inverse
	SEXP WORK = PROTECT(allocVector(REALSXP,dimrtrait));
	F77_CALL(dgetri)(&ntraits, REAL(iZ), &ntraits, INTEGER(IPIV), REAL(WORK), &dimrtrait, &info);
    if (info != 0) {
        if (info > 0) error("The matrix of rates is singular",info);
        error("argument %d had an illegal value",-info);
    }
	// Calculating the matrix vector product
	SEXP S = PROTECT(allocVector(REALSXP,1));
	REAL(S)[0] = 0.0;
	dotprodX(REAL(S), REAL(iZ), REAL(contr), &ntraits, &numbnod);
    
    
    /* Compute the root variance for the Generalized Brownian Motion */
    if(model==4 || model==6 || model==8 || model==9 || model==10){
        
        for(i=0; i<ntraits; i++){
            REAL(muRoot)[i]=(REAL(ancstates)[i]-REAL(mu)[i])/sqrt(REAL(V)[i]);
        }
        /* Dot prod with the root */
        int xroot=1;
        dotprodX(REAL(S), REAL(iZ), REAL(muRoot), &ntraits, &xroot);
    }
    
    
    
	// Calculating sum of the log of contrasts variance for different trees
	REAL(var)[0] = 0.0;
    
        for(i=0; i < nodnbtr; i++){
            REAL(var)[0] += log(REAL(var_contr)[i]);
        }
    
	// fin de var
        for(i=0; i<ntraits; i++){
            REAL(var)[0] += log(REAL(V)[i]);
        }
	
	// resultats
	SEXP result = PROTECT(allocVector(VECSXP, 7));
	SET_VECTOR_ELT(result, 0, contr);
	SET_VECTOR_ELT(result, 1, Z);
	SET_VECTOR_ELT(result, 2, iZ);
	SET_VECTOR_ELT(result, 3, S);
	SET_VECTOR_ELT(result, 4, var);
	SET_VECTOR_ELT(result, 5, det);
	SET_VECTOR_ELT(result, 6, ancstates);

	UNPROTECT(21);
return (result);

}

