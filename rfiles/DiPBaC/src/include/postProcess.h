#ifndef _DiPBaC_POST_PROCESS
#define _DiPBaC_POST_PROCESS

#include <Rcpp.h>

RcppExport SEXP calcDisSimMat(SEXP fileName, SEXP nSweeps, SEXP nBurn, SEXP nFilter,SEXP nSubjects,SEXP nPredictSubjects);

#endif
