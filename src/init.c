#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "LTSM.h"

static R_FortranMethodDef FortRout[] = {
  {"loglik", (DL_FUNC) &F77_SUB(loglik), 44},
  {"dens", (DL_FUNC) &F77_SUB(dens), 44},
  {NULL, NULL, 0}
};


void R_init_LTSM(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
