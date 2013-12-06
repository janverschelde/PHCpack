#ifndef _separate_h
#define _separate_h

#include "DefineType.h"
#include "complex.h"
#include "complexH.h"

// Print functions

void print ( complex<T>* a, complexH<T1>* b, int dim, string var );
void print
 ( complex<T>* a, complexH<T1>** b,
   int dimX, int dimY, int stride, string var );
void print ( complexH<T1>** b, int dimX, int dimY, string var );
void print ( complex<T>* a, int dim, string var );
void print ( complexH<T1>* b, int dim, string var );

// library for the CPU version of the solver

inline complexH<T1> inProd(complexH<T1>* u, complexH<T1>* v, int dim);
inline complexH<T1> inProd_nconj(complexH<T1>* u, complexH<T1>* v, int dim);
void MMprod(complexH<T1>** M1, complexH<T1>** M2, complexH<T1>** PR, int NrM1, int NcM1, int NcM2);
void MVprod(complexH<T1>** M, complexH<T1>* V, complexH<T1>* PR, int NrM, int NcM);
void MVprod_nconj(complexH<T1>** M, complexH<T1>* V, complexH<T1>* PR, int NrM, int NcM);
void matrixTr(complexH<T1>** M, complexH<T1>** Mtr, int dimX, int dimY);
complexH<T1> MatrixError(complexH<T1>** a, complexH<T1>** b, int dimX, int dimY);
complexH<T1> VectorError(complexH<T1>* a, complexH<T1>* b, int dim);
complexH<T1> GPUvsCPUVectorError(complex<T>* a, complexH<T1>* b, int dim);
void CPU_GS_wcopy (complexH<T1>** v, complexH<T1>** v_copy, complexH<T1>** R, int dim, int k);
void CPU_GS(complexH<T1>** v, complexH<T1>** R, int dim, int k);
void matrMult(complexH<T1>** A1, complexH<T1>** A2, int k, int dim);
void BackSubsSec( int dim, complexH<T1>** U, complexH<T1>* y, complexH<T1>* x);

#endif
