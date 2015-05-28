//
//  time_mvmorph.c
//  
//
//  Created by Julien Clavel on 14/11/2014.
//
//

#include "mvmorph.h"



SEXP times_root(SEXP brlength, SEXP edge1, SEXP edge2, SEXP ntip, SEXP Nnode){

    
    int i, nt, ind, e1, e2, nod, ntot;
    
    nt=INTEGER(ntip)[0];
    nod=INTEGER(Nnode)[0];
    ind=nt*2-2;
    ntot=nt+nod;
    
    
    // Edge and alloc vector
    // !! edge must be in postorder or prunningwise order
    PROTECT(edge1 = coerceVector(edge1,INTSXP));
    PROTECT(edge2 = coerceVector(edge2,INTSXP));
    PROTECT(brlength = coerceVector(brlength,REALSXP));
    SEXP times = PROTECT(allocVector(REALSXP,ntot));
    memset(REAL(times),0,(ntot)*sizeof(double));
    
    for(i=ind; i-->0;){
        e2=INTEGER(edge2)[i]-1;
        e1=INTEGER(edge1)[i]-1;
        REAL(times)[e2]=REAL(times)[e1]+REAL(brlength)[i];
    }

    UNPROTECT(4);
    return times;
    
}
