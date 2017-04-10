#include <omp.h>

#include "Southwell.h"

void OMP_MatVecProd(SMEM_MatrixInfo A, double *x, double **y)
{
   
}

void OMP_Residual(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **r)
{
   int n = Mat.A.n;
   #pragma omp parallel num_threads(threads)
   {
      int tid = omp_get_thread_num();
      int nt = Mat.P.part[tid];
      int t_low = Mat.P.disp[tid];
      int k, row, ind;
      double elem; 
      #pragma omp for
      for (int i = 0; i < n; i++){
         (*r)[i] = b[i];
      }
      #pragma omp for
      for (int i = 0; i < n; i++){
         for (int j = 0; j < Mat.A.j_ptr[i+1]-Mat.A.j_ptr[i]; j++){
            ind = Mat.A.j_ptr[i]+j;
            row = Mat.A.i[ind];
            elem = Mat.A.a[ind];
            #pragma omp atomic
            (*r)[row] -= elem * x[i];
         }
      }
   }
}

void OMP_Residual_thread(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **r)
{
   int n = Mat.A.n;
   int tid = omp_get_thread_num();
   int nt = Mat.P.part[tid];
   int t_low = Mat.P.disp[tid];
   int k, row, ind;
   double elem;
   for (int i = 0; i < nt; i++){
      k = t_low + i;
      (*r)[k] = b[k];
   }
   #pragma omp barrier
   for (int i = 0; i < nt; i++){
      k = t_low + i;
      for (int j = 0; j < Mat.A.j_ptr[k+1]-Mat.A.j_ptr[k]; j++){
         ind = Mat.A.j_ptr[k]+j;
         row = Mat.A.i[ind];
         elem = Mat.A.a[ind];
         #pragma omp atomic
         (*r)[row] -= elem * x[k];
      }
   }
}


double OMP_Norm2(
   double *x, 
   int n)
{
   double sum = 0;
   #pragma omp parallel for num_threads(threads) reduction(+:sum)
   for (int i = 0; i < n; i++){
      sum += pow(fabs(x[i]), 2);
   }
   #pragma omp barrier
   return sqrt(sum);
}

void OMP_FreeWriteThread(Write *W)
{
   _mm_free(W->map);
   _mm_free(W->data);
}
