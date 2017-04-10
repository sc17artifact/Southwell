#include "Southwell.h"

void ReorderTriplet_csc(Triplet T, CSC *A, OrderInfo *P)
{
   int row, col;
   int s, q, k, map_i, p;
   double elem;

   int n = A->n;
   int nnz = A->nnz;
   int **col_list = (int **)malloc(n * sizeof(int *));
   double **elem_list = (double **)malloc(n * sizeof(double *));
   int *len = (int *)calloc(n, sizeof(int));
   int *col_flag = (int *)calloc(n, sizeof(int));
   int *p_count = (int *)calloc(P->nparts, sizeof(int));
   A->diag = (double *)calloc(n, sizeof(double));
   A->j_ptr = (int *)calloc(n+1, sizeof(int));
   P->map = (int *)calloc(n, sizeof(int));
   for (int i = 0; i < nnz; i++){
      col = T.j[i];
      row = T.i[i];
      elem = T.a[i];
      p = P->perm[col];
      k = p_count[p];
      if (!col_flag[col]){
         P->map[col] = k + P->disp[p];
         col_flag[col] = 1;
         p_count[p]++;
      }
      map_i = P->map[col];
      len[map_i]++;
      if(row == col) {
         A->diag[map_i] = elem;
      }
   }
   for (int i = 0; i < n; i++){
      col_list[i] = (int *)calloc(len[i], sizeof(int));
      elem_list[i] = (double *)calloc(len[i], sizeof(double));
      len[i] = 0;
   }
   for (int i = 0; i < nnz; i++){
      col = T.j[i];
      row = T.i[i];
      elem = T.a[i];
      p = P->perm[col];
      map_i = P->map[col];
      q = len[map_i];
      col_list[map_i][q] = row;
      elem_list[map_i][q] = elem;
      len[map_i]++;
   }
   A->i = (int *)calloc(nnz,  sizeof(int));
   A->a = (double *)calloc(nnz,  sizeof(double));
   k = 0;
   A->j_ptr[0] = 0;
   for (int i = 0; i < n; i++){
       A->j_ptr[i+1] = A->j_ptr[i] + len[i];
       for (int j = 0; j < len[i]; j++){
          A->i[k] = col_list[i][j];
          A->a[k] = elem_list[i][j];
          k++;
       }
   }
   for (int i = 0; i < n; i++){
      for (int j = 0; j < A->j_ptr[i+1] - A->j_ptr[i]; j++){
         k = A->j_ptr[i] + j;
         A->i[k] = P->map[A->i[k]];
      }
   }
   for (int i = 0; i < n; i++){
      free(col_list[i]);
      free(elem_list[i]);
   }
   free(col_list);
   free(elem_list);
   free(p_count);
   free(col_flag);
   free(len);
}

void Reorder(OrderInfo *P,
             Triplet *T,
             CSC *A)
{
   double start, stop;
   P->disp = (int *)malloc((P->nparts+1) * sizeof(int));
   P->disp[0] = 0;
   for (int i = 0; i < P->nparts; i++){
      P->disp[i+1] = P->disp[i] + P->part[i];
   }
   start = omp_get_wtime();
   ReorderTriplet_csc(*T, A, P);
   stop = omp_get_wtime() - start;
   free(T->a);
   free(T->i);
   free(T->j);
   free(P->perm);
   free(P->map);
   if (!format_out_flag)
      printf("\nMatrix reordered, time = %es.\n", stop);
}

void Metis_to_Triplet(MetisGraph G, Triplet *T)
{
   T->n = (int)G.n;
   T->nnz = (int)G.nnz;
   int q = 0;
   for (int i = 0; i < T->n; i++){
      for (int j = 0; j < (int)(G.xadj[i+1]-G.xadj[i]); j++) {
         T->j[q] = i;
         q++;
      }
   }
   for (int i = 0; i < T->nnz; i++){
      T->i[i] = (int)G.adjncy[i];
      T->a[i] = (double)G.adjwgt[i];
   }
}

void Metis_to_CSC(CSC *A, MetisGraph G)
{
   int q, ind;
   A->n = G.n;
   A->nnz = G.nnz;
   A->diag = (double *)calloc(A->n, sizeof(double));
   A->j_ptr = (int *)malloc((A->n + 1) * sizeof(int));
   A->i = (int *)malloc(A->nnz * sizeof(int));
   A->a = (double *)malloc(A->nnz * sizeof(double));
   for (int i = 0; i < A->n; i++){
      A->j_ptr[i] = G.xadj[i];
      q = (int)(G.xadj[i+1]-G.xadj[i]);
      for (int j = 0; j < q; j++) {
         ind = (int)(G.xadj[i]+j);
         if (i == G.adjncy[ind]){
            A->diag[i] = (double)G.adjwgt[ind];
         }
      }
   }
   A->j_ptr[A->n] = A->nnz;
   for (int i = 0; i < A->nnz; i++){
      A->i[i] = (int)G.adjncy[i];
      A->a[i] = (double)G.adjwgt[i];
   }
}

void Get_Loc_CSC(CSC A, OrderInfo P, CSC *A_loc, int parti)
{
   int q = 0, ind, ai, row;
   int nt = P.part[parti];
   int disp = P.disp[parti];
   double elem;
   A_loc->diag = (double *)calloc(nt, sizeof(double));
   A_loc->j_ptr = (int *)malloc((nt + 1) * sizeof(int));
   A_loc->nnz = 0;
   A_loc->n = nt;
   for (int i = 0; i < nt; i++){
       ai = disp + i;
       A_loc->nnz += A.j_ptr[ai+1] - A.j_ptr[ai];
   }
   A_loc->i = (int *)malloc(A_loc->nnz * sizeof(int));
   A_loc->a = (double *)malloc(A_loc->nnz * sizeof(double));
   A_loc->j_ptr[0] = 0;
   for (int i = 0; i < nt; i++){
       ai = disp + i;
       A_loc->j_ptr[i+1] = A.j_ptr[ai+1] - A.j_ptr[disp];
       for (int j = 0; j < A.j_ptr[ai+1] - A.j_ptr[ai]; j++){
          ind = A.j_ptr[ai]+j;
          row = A.i[ind];
          elem = A.a[ind];
          A_loc->i[q] = row;
          A_loc->a[q] = elem;
          if (ai == row) A_loc->diag[i] = elem;
          q++;
       }
   }
}

void List_to_Block(CSC *A,
                   std::vector<std::list<int>> ind_list,
                   std::vector<std::list<double>> elem_list)
{
   int row_count = 0;
   int prev_row = 0;
   int next_row;
   int k;
   int temp_nnz = 0;

   A->j_ptr = (int *)calloc(A->n+1, sizeof(int));
   A->nnz = 0;
   A->j_ptr[0] = 0;
   for (int i = 0; i < A->n; i++){
      A->nnz += ind_list[i].size();
      A->j_ptr[i+1] = A->j_ptr[i] + ind_list[i].size();
      temp_nnz += (A->j_ptr[i+1] - A->j_ptr[i]);
   }

   A->i = (int *)calloc(A->nnz, sizeof(int));
   A->a = (double *)calloc(A->nnz, sizeof(double));
   k = 0;
   for (int i = 0; i < A->n; i++){
      for (int j = 0; j < A->j_ptr[i+1] - A->j_ptr[i]; j++){
         A->i[k] = ind_list[i].front();
         A->a[k] = elem_list[i].front();
         ind_list[i].pop_front();
         elem_list[i].pop_front();
         k++;
      }
   }
}

void ScaleDiag_CSC(CSC *A)
{
   int row, ind;
   for (int i = 0; i < A->n; i++){
      for (int j = 0; j < A->j_ptr[i+1] - A->j_ptr[i]; j++){
         ind = A->j_ptr[i]+j;
         row = A->i[ind];
         A->a[ind] *= (1/sqrt(A->diag[i]) * 1/sqrt(A->diag[row])); 
      }
   }
   for (int i = 0; i < A->n; i++) A->diag[i] = 1;
}
