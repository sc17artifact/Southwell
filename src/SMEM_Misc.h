#include "Southwell.h"

void OMP_Residual(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **r);

void OMP_Residual_thread(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **r);

double OMP_Norm2(
   double *x,
   int n);
