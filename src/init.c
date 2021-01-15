#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 

   Check these declarations against the C/Fortran source code.

*/

/* .Fortran calls */
extern void F77_NAME(astk)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(circ)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(covst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(listafunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(klistafunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pcffunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(stikfunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gspcore)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gspcoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gtecore)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gtecoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmrcore)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmrcoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmtcore)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmtcoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmmrcore)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmmrcoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmmtcore)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kmmtcoreinh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
static const R_FortranMethodDef FortranEntries[] = {

  {"astk",          (DL_FUNC) &F77_NAME(astk),          15},
  {"circ",          (DL_FUNC) &F77_NAME(circ),          12},
  {"covst",         (DL_FUNC) &F77_NAME(covst),         13},
  {"listafunction", (DL_FUNC) &F77_NAME(listafunction), 27},
  {"klistafunction", (DL_FUNC) &F77_NAME(klistafunction), 23},
  {"pcffunction",   (DL_FUNC) &F77_NAME(pcffunction),   23},
  {"stikfunction",  (DL_FUNC) &F77_NAME(stikfunction),  20},
  {"gspcore",  (DL_FUNC) &F77_NAME(gspcore),  9},
  {"gspcoreinh",  (DL_FUNC) &F77_NAME(gspcoreinh),  16},
  {"gtecore",  (DL_FUNC) &F77_NAME(gtecore),  9},
  {"gtecoreinh",  (DL_FUNC) &F77_NAME(gtecoreinh),  16},
  {"kmrcore",  (DL_FUNC) &F77_NAME(kmrcore),  9},
  {"kmrcoreinh",  (DL_FUNC) &F77_NAME(kmrcoreinh),  16},
  {"kmtcore",  (DL_FUNC) &F77_NAME(kmtcore),  8},
  {"kmtcoreinh",  (DL_FUNC) &F77_NAME(kmtcoreinh),  15},
  {"kmmrcore",  (DL_FUNC) &F77_NAME(kmrcore),  9},
  {"kmmrcoreinh",  (DL_FUNC) &F77_NAME(kmrcoreinh),  16},
  {"kmmtcore",  (DL_FUNC) &F77_NAME(kmtcore),  8},
  {"kmmtcoreinh",  (DL_FUNC) &F77_NAME(kmtcoreinh),  15},
  {NULL, NULL, 0}
};

void R_init_stpp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

