#include <omp.h>

#include "Misc.h"
#include "SMEM_Misc.h"
#include "SMEM_Setup.h"
#include "Southwell.h"

int sync_flag;

int OMP_CheckRelax(
   SMEM_MatrixInfo Mat, 
   SMEM_ThreadInfo *T, 
   double *r, 
   int i, 
   int ai)
{
   return IsZeroMaxAbsDouble(T->N.res[i], T->N.len[i]+1);
}

void OMP_Relax(
   SMEM_MatrixInfo Mat, 
   double **r, 
   double ri, 
   double di, 
   int ai_relax)
{
   for (int j = 0; j < Mat.A.j_ptr[ai_relax+1] -
                       Mat.A.j_ptr[ai_relax]; j++){
      #pragma omp atomic
      (*r)[Mat.A.i[Mat.A.j_ptr[ai_relax]+j]] -= 
         ri * di * Mat.A.a[Mat.A.j_ptr[ai_relax]+j];
   }
}

void OMP_RelaxScatter(
   SMEM_MatrixInfo Mat, 
   double ri, 
   double di, 
   Write *W, 
   int i_relax, 
   int ai_relax)
{
   #pragma simd
   for (int j = 0; j < Mat.A.j_ptr[ai_relax+1] -
                       Mat.A.j_ptr[ai_relax]; j++){
      W->data[W->iter] = Mat.A.a[Mat.A.j_ptr[ai_relax]+j] * di * ri;
      W->iter++;
   }
}

void OMP_RelaxGather(
   SMEM_MatrixInfo Mat, 
   double ri, 
   double di, 
   Write *W, 
   int i_relax, 
   int ai_relax)
{
   #pragma simd
   for (int j = 0; j < Mat.A.j_ptr[ai_relax+1] -
                       Mat.A.j_ptr[ai_relax]; j++){
      W->data[W->map_loc[W->iter]] += 
         Mat.A.a[Mat.A.j_ptr[ai_relax]+j] * di * ri;
      W->iter++;
   }
}

void OMP_SyncParSouthwell(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **u,
   SolveParams params,
   SolveData *out)
{
   sync_flag = 1;
   int n = Mat.A.n;
   double *r = (double *)_mm_malloc(n * sizeof(double), 64);
   double *r_temp = (double *)malloc(n * sizeof(double));
   OMP_Residual(Mat, x, b, &r_temp);
   if ((OMP_Norm2(r_temp,n) < params.tol) ||
       (params.sweep_max == 0)){ 
      for (int i = 0; i < threads; i++){
         out->t[i] = 0;
         out->sweep[i] = 0;
         out->relax[i] = 0;
      }
      for (int i = 0; i < n; i++) out->relax_hist[i] = 0;
      return;
   }
   int converge = 0;
   #pragma omp parallel num_threads(threads)
   {
      SMEM_ThreadInfo T;
      OMP_SetupNeighb_thread(Mat, &T);
      int ai, ind, row, map_i, q;
      int i_low, i_high, ai_low, ai_high, ai_len;
      int sweep = 0, relax = 0;
      double elem, r_norm;
      int tid = omp_get_thread_num();
      int t_low = Mat.P.disp[tid];
      int t_high = Mat.P.disp[tid+1];
      int nt = Mat.P.part[tid];
      int relax_flag = 1;

      int *relax_loc = (int *)calloc(n, sizeof(int));
      double *u_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *r_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *d_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         r[ai] = r_temp[ai];
         u_loc[i] = x[ai];
         r_loc[i] = r[ai];
         d_loc[i] = 1/Mat.A.diag[ai];
      }
      #pragma omp barrier
      double start = omp_get_wtime(), end = 0, elapse;
      #pragma omp barrier
      while(1){
         #pragma omp barrier
         if ((OMP_Norm2(r, n) < params.tol) ||
             (sweep >= params.sweep_max)) break;
         for (int i = 0; i < nt; i++) relax_loc[i] = 0;
         #pragma omp barrier
         start = omp_get_wtime();
         for (int i = 0; i < nt; i++) {
            r_loc[i] = r[t_low+i];
            T.N.res[i][0] = r_loc[i];
         }
         #pragma simd
         for (int i = 0; i < nt; i++){
            #pragma simd
            for (int j = 0; j < T.N.len[i]; j++){
               T.N.res[i][j+1] = r[T.N.rows[i][j]];
            }
         }
         #pragma omp barrier
         for (int i = 0; i < nt; i++){
            ai = t_low + i; ai_low = Mat.A.j_ptr[ai]; ai_high = Mat.A.j_ptr[ai+1];
            ai_len = ai_high - ai_low;
            relax_flag = OMP_CheckRelax(Mat, &T, r, i, ai); 
            if (relax_flag){
               relax_loc[i] = 1;
               OMP_Relax(Mat, &r, r_loc[i], d_loc[i], ai);
               u_loc[i] += r_loc[i] * d_loc[i];
               relax++;
            }
         }
         end = omp_get_wtime();
         elapse += end - start;
         sweep++;
      }
      #pragma omp barrier
      out->t[tid] = elapse;
      out->sweep[tid] = sweep;
      out->relax[tid] = relax;
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         (*u)[ai] = u_loc[i];
         out->relax_hist[ai] = relax_loc[i];
      }
      #pragma omp barrier
      free(relax_loc);
      _mm_free(u_loc);
      _mm_free(r_loc);
      _mm_free(d_loc);
   }
   _mm_free(r);
   free(r_temp);
}

void OMP_SyncJacobi(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **u,
   SolveParams params,
   SolveData *out)

{
   sync_flag = 1;
   int n = Mat.A.n;
   double *r = (double *)_mm_malloc(n * sizeof(double), 64);
   double *r_temp = (double *)malloc(n * sizeof(double));
   OMP_Residual(Mat, x, b, &r_temp);
   if ((OMP_Norm2(r_temp,n) < params.tol) ||
       (params.sweep_max == 0)){
      for (int i = 0; i < threads; i++){
         out->t[i] = 0;
         out->sweep[i] = 0;
         out->relax[i] = 0;
      }
      return;
   }
   int converge = 0;
   #pragma omp parallel num_threads(threads)
   {
      Write W;
      OMP_SetupWrite_gather(Mat, &W);
      int ai, ind, row, map_i, q;
      int i_low, i_high, ai_low, ai_high, ai_len;
      int sweep = 0, relax = 0, comm = 0;
      double elem, r_norm;
      int tid = omp_get_thread_num();
      int t_low = Mat.P.disp[tid];
      int t_high = Mat.P.disp[tid+1];
      int nt = Mat.P.part[tid];
      int relax_flag = 1;

      double *u_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *r_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *d_loc = (double *)_mm_malloc(nt * sizeof(double), 64);

      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         r[ai] = r_temp[ai];
         u_loc[i] = x[ai];
         r_loc[i] = r[ai];
         d_loc[i] = 1/Mat.A.diag[ai];
      }
      double start = omp_get_wtime(), end = 0, elapse;
      #pragma omp barrier
      while(1){
         if ((OMP_Norm2(r, n) < params.tol) ||
             (sweep >= params.sweep_max)) break;
         #pragma omp barrier       
         W.iter = 0;
         start = omp_get_wtime();
         for (int i = 0; i < nt; i++){
            ai = t_low + i; ai_low = Mat.A.j_ptr[ai]; ai_high = Mat.A.j_ptr[ai+1];
            ai_len = ai_high - ai_low;
            OMP_RelaxGather(Mat, r_loc[i], d_loc[i], &W, i, ai);
            relax++;
            comm += ai_len;
         }
         for (int i = 0; i < nt; i++) u_loc[i] += r_loc[i] * d_loc[i];
         for (int i = 0; i < W.data_len; i++){
            #pragma omp atomic
            r[W.map[i]] -= W.data[i];
            W.data[i] = 0;
         }
         #pragma omp barrier
         for (int i = 0; i < nt; i++) r_loc[i] = r[t_low+i];
         end = omp_get_wtime();
         elapse += end - start;
         sweep++;
      }
      #pragma omp barrier
      out->t[tid] = elapse;
      out->sweep[tid] = sweep;
      out->relax[tid] = relax;
      out->comm[tid] = comm;
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         (*u)[ai] = u_loc[i];
      }
      #pragma omp barrier
      _mm_free(u_loc);
      _mm_free(r_loc);
      _mm_free(d_loc);
   }
   _mm_free(r);
   free(r_temp);
}


void OMP_AsyncParSouthwell(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **u,
   SolveParams params,
   SolveData *out)

{
   sync_flag = 1;
   int n = Mat.A.n;
   double *r = (double *)_mm_malloc(n * sizeof(double), 64);
   double *r_temp = (double *)malloc(n * sizeof(double));
   int *converge_flag = (int *)calloc(threads, sizeof(int));
   OMP_Residual(Mat, x, b, &r_temp);
   if ((OMP_Norm2(r_temp,n) < params.tol) ||
       (params.sweep_max == 0)){
      for (int i = 0; i < threads; i++){
         out->t[i] = 0;
         out->sweep[i] = 0;
         out->relax[i] = 0;
      }
      return;
   }
   int converge = 0;
   #pragma omp parallel num_threads(threads)
   {
      SMEM_ThreadInfo T;
      OMP_SetupNeighb_thread(Mat, &T);
      int ai, ind, row, map_i;
      int i_low, i_high, ai_low, ai_high, ai_len;
      int sweep = 0, relax = 0, comm = 0, q;
      double elem, r_norm;
      int tid = omp_get_thread_num();
      int t_low = Mat.P.disp[tid];
      int t_high = Mat.P.disp[tid+1];
      int nt = Mat.P.part[tid];
      int relax_flag = 1;

      double *u_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *r_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *d_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         r[ai] = r_temp[ai];
         u_loc[i] = x[ai];
         r_loc[i] = r[ai];
         d_loc[i] = 1/Mat.A.diag[ai];
      }
      double start = omp_get_wtime(), end = 0, elapse;
      #pragma omp barrier
      while(1){
         if (!converge_flag[tid]){
            if ((Norm2(r_loc, nt) < params.tol) ||
                (sweep >= params.sweep_max)){
               converge_flag[tid] = 1;
               #pragma omp critical
               converge++;
            }
         }
         if (converge >= threads) break;
         start = omp_get_wtime();
         for (int i = 0; i < nt; i++) T.N.res[i][0] = r_loc[i];
         #pragma simd
         for (int i = 0; i < nt; i++){
            #pragma simd
            for (int j = 0; j < T.N.len[i]; j++){
               T.N.res[i][j+1] = r[T.N.rows[i][j]];
            }
         }
         for (int i = 0; i < nt; i++){
            ai = t_low + i; ai_low = Mat.A.j_ptr[ai]; ai_high = Mat.A.j_ptr[ai+1];
            ai_len = ai_high - ai_low;
            relax_flag = OMP_CheckRelax(Mat, &T, r, i, ai); 
            if (relax_flag){
               OMP_Relax(Mat, &r, r_loc[i], d_loc[i], ai);
               u_loc[i] += r_loc[i] * d_loc[i];
               relax++;
            }
         }
         for (int i = 0; i < nt; i++) r_loc[i] = r[t_low+i];
         end = omp_get_wtime();
         elapse += end - start;
         sweep++;
      }
      #pragma omp barrier
      out->t[tid] = elapse;
      out->sweep[tid] = sweep;
      out->relax[tid] = relax;
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         (*u)[ai] = u_loc[i];
      }
      #pragma omp barrier
      _mm_free(u_loc);
      _mm_free(r_loc);
      _mm_free(d_loc);
   }
   _mm_free(r);
   free(r_temp);
   free(converge_flag);
}


void OMP_AsyncJacobi(
   SMEM_MatrixInfo Mat, 
   double *x, 
   double *b, 
   double **u,
   SolveParams params,
   SolveData *out)
{
   sync_flag = 1;
   int n = Mat.A.n;
   double *r = (double *)_mm_malloc(n * sizeof(double), 64);
   double *r_temp = (double *)malloc(n * sizeof(double));
   int *converge_flag = (int *)calloc(threads, sizeof(int));
   OMP_Residual(Mat, x, b, &r_temp);
   if ((OMP_Norm2(r_temp,n) < params.tol) ||
       (params.sweep_max == 0)){
      for (int i = 0; i < threads; i++){
         out->t[i] = 0;
         out->sweep[i] = 0;
         out->relax[i] = 0;
      }
      return;
   }
   int converge = 0;
   #pragma omp parallel num_threads(threads)
   {
      Write W;
      OMP_SetupWrite_gather(Mat, &W);
      int ai, ind, row, map_i;
      int i_low, i_high, ai_low, ai_high, ai_len;
      int sweep = 0, relax = 0, comm = 0;
      double elem, r_norm;
      int tid = omp_get_thread_num();
      int t_low = Mat.P.disp[tid];
      int t_high = Mat.P.disp[tid+1];
      int nt = Mat.P.part[tid];
      int relax_flag = 1;

      double *u_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *r_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      double *d_loc = (double *)_mm_malloc(nt * sizeof(double), 64);
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         r[ai] = r_temp[ai];
         u_loc[i] = x[ai];
         r_loc[i] = r[ai];
         d_loc[i] = 1/Mat.A.diag[ai];
      }
      double start = omp_get_wtime(), end = 0, elapse;
      #pragma omp barrier
      while(1){
         if (!converge_flag[tid]){
            if ((Norm2(r_loc, nt) < params.tol) ||
                (sweep >= params.sweep_max)){
               converge_flag[tid] = 1;
               #pragma omp critical
               converge++;
            }
         }
         if (converge >= threads) break;
         W.iter = 0;
         start = omp_get_wtime();
         for (int i = 0; i < nt; i++){
            ai = t_low + i; ai_low = Mat.A.j_ptr[ai]; ai_high = Mat.A.j_ptr[ai+1];
            ai_len = ai_high - ai_low;
            OMP_RelaxGather(Mat, r_loc[i], d_loc[i], &W, i, ai);
            relax++;
            comm += ai_len;
         }
         for (int i = 0; i < nt; i++) u_loc[i] += r_loc[i] * d_loc[i];
         for (int i = 0; i < W.data_len; i++){
            #pragma omp atomic
            r[W.map[i]] -= W.data[i];
            W.data[i] = 0;
         }
         for (int i = 0; i < nt; i++) r_loc[i] = r[t_low+i];
         end = omp_get_wtime();
         elapse += end - start;
         sweep++;
      }
      #pragma omp barrier
      out->t[tid] = elapse;
      out->sweep[tid] = sweep;
      out->relax[tid] = relax;
      out->comm[tid] = comm;
      for (int i = 0; i < nt; i++){
         ai = t_low + i;
         (*u)[ai] = u_loc[i];
      }
      #pragma omp barrier
      _mm_free(u_loc);
      _mm_free(r_loc);
      _mm_free(d_loc);
   }
   _mm_free(r);
   free(r_temp);
   free(converge_flag);
}

void OMP_MCGS(
   SMEM_MatrixInfo Mat,
   double *x,
   double *b, 
   double **u,
   SolveParams params,
   SolveData *out)

{
   sync_flag = 0;
   int n = Mat.A.n;
   double *r = (double *)calloc(n, sizeof(double));
   double *r_temp = (double *)calloc(n, sizeof(double));
   OMP_Residual(Mat, x, b, &r_temp);
   if ((OMP_Norm2(r_temp,n) < params.tol) ||
       (params.sweep_max == 0)){
      for (int i = 0; i < threads; i++){
         out->t[i] = 0;
         out->sweep[i] = 0;
         out->relax[i] = 0;
      }
      return;
   }
   #pragma omp parallel num_threads(threads)
   {
      int ind, row;
      double elem, rj, dj;
      int sweep = 0, relax = 0, comm = 0;
      int tid = omp_get_thread_num();
      int relax_flag = 1;
      #pragma omp for schedule(guided)
      for (int i = 0; i < n; i++){
         (*u)[i] = x[i];
         r[i] = r_temp[i];
      }
      double start = omp_get_wtime(), end = 0, elapse;
      #pragma omp barrier
      while(1){
         if ((OMP_Norm2(r, n) < params.tol) ||
             (sweep >= params.sweep_max)) break; 
         #pragma omp barrier
         for (int i = 0; i < Mat.P.nparts; i++){
            start = omp_get_wtime();
            #pragma omp for schedule (guided)
            for (int j = Mat.P.disp[i]; j < Mat.P.disp[i]+Mat.P.part[i]; j++){
               rj = r[j];
               dj = 1. / Mat.A.diag[j];
               (*u)[j] += rj * dj;
               OMP_Relax(Mat, &r, rj, dj, j);
               relax++;
            }
            end = omp_get_wtime();
            elapse += end - start;
            sweep++;
            if (sweep > params.sweep_max) break;
         }
      }
      #pragma omp barrier
      out->t[tid] = elapse;
      out->sweep[tid] = sweep;
      out->relax[tid] = relax;
   }
   free(r);
   free(r_temp);
}
