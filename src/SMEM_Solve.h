#include "Southwell.h"

void OMP_SyncParSouthwell(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **u,
   SolveParams params,
   SolveData *out);

void OMP_SyncJacobi(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **u,
   SolveParams params,
   SolveData *out);

void OMP_AsyncParSouthwell(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **u,
   SolveParams params,
   SolveData *out);

void OMP_AsyncJacobi(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **u,
   SolveParams params,
   SolveData *out);

void OMP_MCGS(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b,
   double **u,
   SolveParams params,
   SolveData *out);
