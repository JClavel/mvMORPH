#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP Chol_RPF(SEXP A, SEXP D, SEXP dat, SEXP nterm, SEXP ndimA, SEXP mserr, SEXP ismserr);
extern SEXP Chol_RPF_only(SEXP A, SEXP ndimA, SEXP mserr, SEXP ismserr);
extern SEXP Chol_RPF_quadprod(SEXP U, SEXP resid, SEXP nterm);
extern SEXP Chol_RPF_quadprod_column(SEXP U, SEXP resid, SEXP nterm);
extern SEXP Chol_RPF_univ(SEXP A, SEXP D, SEXP dat, SEXP nterm, SEXP ndimA, SEXP mserr, SEXP ismserr);
extern SEXP Chol_RPF_univ_only(SEXP A, SEXP ndimA, SEXP mserr, SEXP ismserr);
extern SEXP Expect_matrix(SEXP S1, SEXP S, SEXP lambda, SEXP time, SEXP theta0, SEXP theta1, SEXP matdiag);
extern SEXP givens_ortho (SEXP Q, SEXP angle, SEXP ndim);
extern SEXP kronecker_mvmorph(SEXP R, SEXP C, SEXP Rrows, SEXP Crows);
extern SEXP kronecker_shift(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP V);
extern SEXP kronecker_shiftEB_BM(SEXP R1, SEXP R2, SEXP C1, SEXP C2, SEXP beta, SEXP Rrows, SEXP Crows);
extern SEXP kronecker_shiftEB_OU(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP V);
extern SEXP kroneckerEB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows);
extern SEXP kroneckerSpar_shift(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP V, SEXP IA, SEXP JA, SEXP A);
extern SEXP kroneckerSpar_shift_EB_BM(SEXP R1, SEXP R2, SEXP C1, SEXP C2, SEXP beta, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A);
extern SEXP kroneckerSpar_shift_OU_EB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP V, SEXP IA, SEXP JA, SEXP A);
extern SEXP kroneckerSparEB(SEXP R, SEXP C, SEXP beta, SEXP Rrows, SEXP Crows, SEXP IA, SEXP JA, SEXP A);
extern SEXP kroneckerSum(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP dimlist);
extern SEXP kroneckerSumSpar(SEXP R, SEXP C, SEXP Rrows, SEXP Crows, SEXP dimlist, SEXP IA, SEXP JA, SEXP A);
extern SEXP mvmorph_covar_mat (SEXP nterm, SEXP bt,SEXP lambda, SEXP S, SEXP sigmasq, SEXP S1);
extern SEXP mvmorph_covar_ou_fixed(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP mvmorph_covar_ou_random(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP mvmorph_covar_ou_rpf_fixed(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP mvmorph_covar_ou_rpf_random(SEXP A, SEXP alpha, SEXP sigma);
extern SEXP mvmorph_covar_ou_sparse (SEXP A, SEXP JA, SEXP IA, SEXP nterm, SEXP bt,SEXP lambda, SEXP S, SEXP sigmasq, SEXP S1);
extern SEXP mvmorph_weights (SEXP nterm, SEXP epochs, SEXP lambda, SEXP S, SEXP S1, SEXP beta, SEXP root);
extern SEXP PIC_gen(SEXP x, SEXP n, SEXP Nnode, SEXP nsp, SEXP edge1, SEXP edge2, SEXP edgelength, SEXP times, SEXP rate, SEXP Tmax, SEXP Model, SEXP mu, SEXP sigma);
extern SEXP seq_root2tipM(SEXP edge, SEXP nbtip, SEXP nbnode);
extern SEXP simmap_covar (SEXP nterm, SEXP bt, SEXP lambda, SEXP S, SEXP S1, SEXP sigmasq);
extern SEXP spherical(SEXP param, SEXP variance, SEXP dim);
extern SEXP squareRootM(SEXP edge1, SEXP edge2, SEXP edgelength, SEXP nsp, SEXP inverse);
extern SEXP times_root(SEXP brlength, SEXP edge1, SEXP edge2, SEXP ntip, SEXP Nnode);
extern SEXP Weight_matrix(SEXP S1, SEXP S, SEXP lambda, SEXP time, SEXP matdiag);

static const R_CallMethodDef CallEntries[] = {
    {"Chol_RPF",                    (DL_FUNC) &Chol_RPF,                     7},
    {"Chol_RPF_only",               (DL_FUNC) &Chol_RPF_only,                4},
    {"Chol_RPF_quadprod",           (DL_FUNC) &Chol_RPF_quadprod,            3},
    {"Chol_RPF_quadprod_column",    (DL_FUNC) &Chol_RPF_quadprod_column,     3},
    {"Chol_RPF_univ",               (DL_FUNC) &Chol_RPF_univ,                7},
    {"Chol_RPF_univ_only",          (DL_FUNC) &Chol_RPF_univ_only,           4},
    {"Expect_matrix",               (DL_FUNC) &Expect_matrix,                7},
    {"givens_ortho",                (DL_FUNC) &givens_ortho,                 3},
    {"kronecker_mvmorph",           (DL_FUNC) &kronecker_mvmorph,            4},
    {"kronecker_shift",             (DL_FUNC) &kronecker_shift,              5},
    {"kronecker_shiftEB_BM",        (DL_FUNC) &kronecker_shiftEB_BM,         7},
    {"kronecker_shiftEB_OU",        (DL_FUNC) &kronecker_shiftEB_OU,         6},
    {"kroneckerEB",                 (DL_FUNC) &kroneckerEB,                  5},
    {"kroneckerSpar_shift",         (DL_FUNC) &kroneckerSpar_shift,          8},
    {"kroneckerSpar_shift_EB_BM",   (DL_FUNC) &kroneckerSpar_shift_EB_BM,   10},
    {"kroneckerSpar_shift_OU_EB",   (DL_FUNC) &kroneckerSpar_shift_OU_EB,    9},
    {"kroneckerSparEB",             (DL_FUNC) &kroneckerSparEB,              8},
    {"kroneckerSum",                (DL_FUNC) &kroneckerSum,                 5},
    {"kroneckerSumSpar",            (DL_FUNC) &kroneckerSumSpar,             8},
    {"mvmorph_covar_mat",           (DL_FUNC) &mvmorph_covar_mat,            6},
    {"mvmorph_covar_ou_fixed",      (DL_FUNC) &mvmorph_covar_ou_fixed,       3},
    {"mvmorph_covar_ou_random",     (DL_FUNC) &mvmorph_covar_ou_random,      3},
    {"mvmorph_covar_ou_rpf_fixed",  (DL_FUNC) &mvmorph_covar_ou_rpf_fixed,   3},
    {"mvmorph_covar_ou_rpf_random", (DL_FUNC) &mvmorph_covar_ou_rpf_random,  3},
    {"mvmorph_covar_ou_sparse",     (DL_FUNC) &mvmorph_covar_ou_sparse,      9},
    {"mvmorph_weights",             (DL_FUNC) &mvmorph_weights,              7},
    {"PIC_gen",                     (DL_FUNC) &PIC_gen,                     13},
    {"seq_root2tipM",               (DL_FUNC) &seq_root2tipM,                3},
    {"simmap_covar",                (DL_FUNC) &simmap_covar,                 6},
    {"spherical",                   (DL_FUNC) &spherical,                    3},
    {"squareRootM",                 (DL_FUNC) &squareRootM,                  5},
    {"times_root",                  (DL_FUNC) &times_root,                   5},
    {"Weight_matrix",               (DL_FUNC) &Weight_matrix,                5},
    {NULL, NULL, 0}
};

void R_init_mvMORPH(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
