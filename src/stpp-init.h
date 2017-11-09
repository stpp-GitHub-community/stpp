#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(circ)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(covst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(listafunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pcffunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(stikfunction)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"circ",         (DL_FUNC) &F77_NAME(circ),         12},
  {"covst",        (DL_FUNC) &F77_NAME(covst),        13},
  {"listafunction", (DL_FUNC) &F77_NAME(listafunction), 27},
  {"pcffunction",  (DL_FUNC) &F77_NAME(pcffunction),  23},
  {"stikfunction", (DL_FUNC) &F77_NAME(stikfunction), 20},
  {NULL, NULL, 0}
};

void R_init_stpp(DllInfo *dll);