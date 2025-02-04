#include <omp.h>

#include "Southwell.h"
#include "FileUtils.h"
#include "MatrixUtils.h"
#include "Misc.h"
#include "Multicolor.h"
#include "Laplace.h"

void SMEM_ReorderTriplet_csc(Triplet T, SMEM_MatrixInfo *Mat);
void SMEM_Multicolor(SMEM_MatrixInfo *Mat);

void SMEM_Setup(
   FILE *in_file, 
   SMEM_MatrixInfo *Mat, 
   SMEM_ThreadInfo *Threads,
   int m, 
   int w)
{
   idx_t ncon = 1, objval, n;
   int flag = METIS_OK, rank_p = 0, num_p = 1, q, ind, nnz;
   int imbalance = 0;
   MetisGraph G;
   Triplet T;
   idx_t nparts;
   if (color_flag){
      nparts = 2;
   }
   else {
      nparts = (idx_t)(threads) + imbalance;
   }
   
   idx_t options[METIS_NOPTIONS];
   char buffer[256];
   FILE *temp_file;
   struct timeval tv_start, tv_stop, tv_diff;

   if ((rank_p == 0) && (!format_out_flag))
      printf("\n**************** INITIALIZING PROBLEM ***************\n");
   gettimeofday(&tv_start, NULL);
   if (mat_file_flag){
      ReadBinary_fread_metis(in_file, &G, &T);
   }
   else {
      Laplace2D_FD5pt_metis(&G, &T, m, w);
   }

   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

   n = G.n;
   nnz = (int)G.two_m;
   Mat->A.n = (int)n;
   Mat->A.nnz = nnz;
   idx_t *perm = (idx_t *)calloc((int)G.n, sizeof(idx_t));
   Mat->P.perm = (int *)malloc(n * sizeof(int));
   gettimeofday(&tv_stop, NULL);
   timersub(&tv_stop, &tv_start, &tv_diff);
   if ((rank_p == 0) && (!format_out_flag)){
      printf("\nMATRIX LOADED, time = %ld.%06ld\n",
             (long int)tv_diff.tv_sec, (long int)tv_diff.tv_usec);
      printf("Test problem: ");
      if (mat_file_flag){
         printf("reading matrix from file\n");
      }
      else {
         printf("finite difference laplace\n");
      }
      printf("nnz = %d, n = %d\n", Mat->A.nnz, Mat->A.n);
   }
   if (nparts > 1){
      if (color_flag){
         Metis_to_CSC(&Mat->A, G);
         if (!format_out_flag) printf("\nREORDERING USING MULTICOLORING\n");
         Multicolor(Mat->A, &(Mat->P));
      }
      else { 
         if (nparts > 1){
         if (!format_out_flag) printf("\nREORDERING USING METIS\n");
         flag =  METIS_PartGraphKway(&n, &ncon, G.xadj, G.adjncy, NULL, NULL,
                                     NULL, &nparts, NULL, NULL, options,
                                     &objval, perm);
         }
         if (flag != METIS_OK)
            printf("****WARNING****: METIS returned error\n");
 
         Mat->P.nparts = threads; 
         Mat->P.part = (int *)calloc(Mat->P.nparts, sizeof(int));
         for (int i = 0; i < (int)n; i++){
            Mat->P.perm[i] = (int)perm[i];
            Mat->P.part[Mat->P.perm[i]]++;
            /* IMBALANCE */
           // if ((int)perm[i] > threads-1){
           //    Mat->P.perm[i] = threads-1;
           //    Mat->P.part[Mat->P.perm[i]]++;
           // }
           // else {
           //    Mat->P.perm[i] = (int)perm[i];
           //    Mat->P.part[Mat->P.perm[i]]++;
           // }
         }
      }
      gettimeofday(&tv_stop, NULL);
      timersub(&tv_stop, &tv_start, &tv_diff);
      if (!format_out_flag) printf("GRAPH PARTITIONED, time = %ld.%06ld\n",
             (long int)tv_diff.tv_sec, (long int)tv_diff.tv_usec);
      if (color_flag){
         if (!format_out_flag){ 
             printf("colors = %d\n" 
                    "partition: mean = %.2f, max %d, min %d\n",
                    Mat->P.nparts, MeanInt(Mat->P.part, Mat->P.nparts), 
                    MaxInt(Mat->P.part,nparts), MinInt(Mat->P.part,nparts));
             printf("%d & %d & %d & %d & %d\\\\", 
                    nnz, n, Mat->P.nparts, MaxInt(Mat->P.part,nparts), 
                    MinInt(Mat->P.part,nparts));
         }
      }
      Mat->P.disp = (int *)malloc((Mat->P.nparts+1) * sizeof(int));
      Mat->P.disp[0] = 0;
      for (int i = 1; i <= Mat->P.nparts; i++){
         Mat->P.disp[i] = Mat->P.disp[i-1] + Mat->P.part[i-1];
      }
      gettimeofday(&tv_start, NULL);
      if (reorder_flag){
         ReorderTriplet_csc(T, &(Mat->A), &(Mat->P)); 
         free(T.a);
         free(T.i);
         free(T.j);
         gettimeofday(&tv_stop, NULL);
         timersub(&tv_stop, &tv_start, &tv_diff);
         if (rank_p == 0){
             if (!format_out_flag) printf("\nMATRIX REORDERED, time = %ld.%06ld\n",
                    (long int)tv_diff.tv_sec, (long int)tv_diff.tv_usec);
         }
      }
      else Metis_to_CSC(&Mat->A, G);
   }
   else{
      Metis_to_CSC(&Mat->A, G);
   }
   Write_csc(Mat->A, 1);
   if (!format_out_flag) printf("\n******************************************************\n\n");
   free(perm);
   free(Mat->P.perm);
}


//void OMP_GetBoundary(CSC_Block B,  OMP_Bound *bound)
//{
//   int n = B.n;
//   int nnz = B.nnz;
//   int row;
//   int *row_flag = (int *)calloc(n, sizeof(int));
//   bound->len = 0;
//   for (int i = 0; i < nnz; i++){
//      row = B.i[i];
//      if (!row_flag[row]){
//         row_flag[row] = 1;
//         bound->len++;
//      }
//   }
//   bound->row = (int *)calloc(bound->len, sizeof(int));
//   int q = 0;
//   for (int i = 0; i < n; i++){
//      if (row_flag[i]){
//         bound->row[q] = i;
//         q++;
//      }
//   }
//   free(row_flag);
//}


void OMP_SetupNeighb_thread(SMEM_MatrixInfo Mat, SMEM_ThreadInfo *T)
{
   int tid = omp_get_thread_num();
   SetupNeighb(Mat.A, Mat.P, &T->N, tid);
}

void OMP_SetupWrite_scatter(SMEM_MatrixInfo Mat, Write *W)
{
   int temp_len = 0;
   int row, ind, ai;
   int tid = omp_get_thread_num();
   int nt = Mat.P.part[tid];
   int t_low = Mat.P.disp[tid];
   W->data_len = 0;
   for (int i = 0; i < nt; i++){
      ai = t_low + i;
      for (int j = 0; j < Mat.A.j_ptr[ai+1] - Mat.A.j_ptr[ai]; j++) W->data_len++;
   }
   W->map_len = W->data_len;
   W->data = (double *)_mm_malloc(W->data_len * sizeof(double), 64);
   W->map = (int *)_mm_malloc(W->data_len * sizeof(int), 64);
   W->iter = 0;
   for (int i = 0; i < nt; i++){
      ai = t_low + i;
      for (int j = 0; j < Mat.A.j_ptr[ai+1]-Mat.A.j_ptr[ai]; j++) {
         ind = Mat.A.j_ptr[ai]+j;
         row = Mat.A.i[ind];
         W->map[W->iter] = row;
         W->iter++;
      }
   }
}

void OMP_SetupWrite_gather(SMEM_MatrixInfo Mat, Write *W)
{
   int temp_len = 0;
   int row, ind, ai;
   int tid = omp_get_thread_num();
   int nt = Mat.P.part[tid];
   int t_low = Mat.P.disp[tid];
   int q;
   int n = Mat.A.n;
   W->low = n;
   W->high = 0;
   W->map_len = 0;
   W->data_len = 0;

   int *row_flag = (int *)calloc(n, sizeof(int));

   for (int i = 0; i < nt; i++){
      ai = t_low + i;
      for (int j = 0; j < Mat.A.j_ptr[ai+1]-Mat.A.j_ptr[ai]; j++) {
         ind = Mat.A.j_ptr[ai]+j;
         row = Mat.A.i[ind];
         W->map_len++;
         if (!row_flag[row]){
            W->data_len++;
            row_flag[row] = 1;
         }
      }
   }
   W->data = (double *)_mm_malloc(W->data_len*sizeof(double), 64);
   W->map = (int *)_mm_malloc(W->data_len*sizeof(int), 64);
   W->map_loc = (int *)_mm_malloc(W->map_len*sizeof(int), 64);
   W->iter = 0;
   for (int i = 0; i < n; i++){
      if (row_flag[i]){
         W->map[W->iter] = i;
         W->data[W->iter] = 0;
         row_flag[i] = W->iter;
         W->iter++;
      }
   }
   W->iter = 0;
   for (int i = 0; i < nt; i++){
      ai = t_low + i;
      for (int j = 0; j < Mat.A.j_ptr[ai+1]-Mat.A.j_ptr[ai]; j++) {
         ind = Mat.A.j_ptr[ai]+j;
         row = Mat.A.i[ind];
         W->map_loc[W->iter] = row_flag[row];
         W->iter++;
      }
   }
   free(row_flag);
}
