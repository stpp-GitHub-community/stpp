#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */

extern void F77_NAME(astk)(double *x, double *y, double *txy, int *n, double *lambda, double *ag, double *s, int *ns, double *t, int *nt, double *wbi, double *wbimod, double *wt, int *correc, double *astkf);
extern void F77_NAME(circ)(double *CEM, int *M, int *N, double *xlim, double *ylim, double *tlim, double *model, double *param, double *sigma2, double *scale, double *aniso, double *ani);
extern void F77_NAME(covst)(double *gs, double *xx, double *yy, double *tt, int *nx, int *ny, int *nt, int *model, double *param, double *sigma2, double *scale, double *aniso, double *ani);
extern void F77_NAME(listafunction)(int *i, double *xi, double *yi, double *ti, double *x, double *y, double *txy, int *n, double *xp, double *yp, int *np, double *s, int *ns, double *t, int *nt, double *bsupt, double *binft, double *lambda, int *ks, int *kt, double *hs, double *ht, double *listahat, double *wbi, double *wbimod, double *wt, int *correc);
extern void F77_NAME(pcffunction)(double *x, double *y, double *txy, int *n, double *xp, double *yp, int *np, double *s, int *ns, double *t, int *nt, double *bsupt, double *binft, double *lambda, int *ks, int *kt, double *hs, double *ht, double *pcfhat, double *wbi, double *wbimod, double *wt, int *correc);
extern void F77_NAME(stikfunction)(double *x, double *y, double *txy, int *n, double *xp, double *yp, int *np, double *s, int *ns, double *t, int *nt, double *bsupt, double *binft, double *lambda, int *infd, double *hkhat, double *wbi, double *wbimod, double *wt, int *correc);

static const R_FortranMethodDef FortranEntries[] = {
  {"astk",          (DL_FUNC) &F77_NAME(astk),          15},
  {"circ",          (DL_FUNC) &F77_NAME(circ),          12},
  {"covst",         (DL_FUNC) &F77_NAME(covst),         13},
  {"listafunction", (DL_FUNC) &F77_NAME(listafunction), 27},
  {"pcffunction",   (DL_FUNC) &F77_NAME(pcffunction),   23},
  {"stikfunction",  (DL_FUNC) &F77_NAME(stikfunction),  20},
  {NULL, NULL, 0}
};

void R_init_stpp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
