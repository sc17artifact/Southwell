#include "Southwell.h"
#include "MatrixUtils.h"

void Laplace2D_FD5pt_metis(
   MetisGraph *G, 
   Triplet *T,
   int m, 
   int n)
{
   int col;
   int N = m*n;

   G->n = N;
   G->nnz = (idx_t)(m * (n + 2*(n - 1)) + 2*(m-1) * n);
   G->xadj = (idx_t *)calloc(N, sizeof(idx_t));
   G->adjncy = (idx_t *)calloc((int)G->nnz, sizeof(idx_t));
   G->adjwgt = (real_t *)calloc((int)G->nnz, sizeof(real_t));

   T->nnz = (int)G->nnz;
   T->n = N;
   T->i = (int *)calloc(T->nnz, sizeof(int));
   T->j = (int *)calloc(T->nnz, sizeof(int));
   T->a = (double *)calloc(T->nnz, sizeof(double));
  
   int block_end = m-1;
   int block_start = 0;
   int k = 0;
   G->xadj[0] = 0;
   for(int i = 0; i < N; i++){
      col = i - m;
      if (col >= 0){
         G->adjwgt[k] = -1;
         G->adjncy[k] = col;
         k++;
      }
      if (i > block_start){
         col = i - 1;
         G->adjwgt[k] = -1;
         G->adjncy[k] = col;
         k++;
      }
      G->adjwgt[k] = 4;
      G->adjncy[k] = i;
      k++;
      if (i < block_end){
         col = i + 1;
         G->adjwgt[k] = -1;
         G->adjncy[k] = col;
         k++;
      }
      col = i + m;
      if (col < N){
         G->adjwgt[k] = -1;
         G->adjncy[k] = col;
         k++;
      }
      G->xadj[i+1] = k;

      if (i == block_end){
         block_end += m;
         block_start += m;
      }
   }
   Metis_to_Triplet(*G, T);
}
