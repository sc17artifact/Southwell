#include "Southwell.h"
#include "FileUtils.h"
#include "MatrixUtils.h"
#include "Misc.h"
#include "Multicolor.h"
#include "Laplace.h"

using namespace std;

extern int block_flag;
extern int loc_direct_flag;

void SEQ_SetupBlocks(CSC A,
                     OrderInfo P,
                     CSC **D,
                     CSC **B);

void SEQ_Setup(FILE *in_file, 
               CSC *A, 
               CSC **B,
               CSC **D,
               OrderInfo *P,
               PardisoInfo **Pard,
               int nx, 
               int ny)
{
   idx_t n, objval, ncon = 1;
   int nnz;
   int flag;
   MetisGraph G;
   Triplet T;
   double start, stop;

   char buffer[256];
   FILE *temp_file;

   if (!format_out_flag)
      printf("\n**************** INITIALIZING PROBLEM ***************\n");
   start = omp_get_wtime();
   if (mat_file_flag){
      ReadBinary_fread_metis(in_file, &G, &T);
   }
   else {
      Laplace2D_FD5pt_metis(&G, &T, nx, ny);
   }

   n = G.n;
   nnz = (int)G.nnz;
   A->n = (int)n;
   A->nnz = nnz;
   P->perm = (int *)malloc(n * sizeof(int));
   if (!format_out_flag){
      printf("\nMATRIX LOADED, time = %e\n",
             omp_get_wtime() - start);
      printf("Test problem: ");
      if (mat_file_flag){
         printf("reading matrix from file\n");
      }
      else {
         printf("finite difference laplace\n");
      }
      printf("& nnz & n &\n");      
      printf("& %d & %d &\n", A->nnz, A->n);
   }
   Metis_to_CSC(A, G);
   ScaleDiag_CSC(A);
   if (color_flag){
      start = omp_get_wtime();
      if (!format_out_flag) printf("\nREORDERING USING MULTICOLORING\n");
      Multicolor(*A, P);
      P->disp = (int *)malloc((P->nparts+1) * sizeof(int));
      P->disp[0] = 0;
      for (int i = 0; i < P->nparts; i++){
         P->disp[i+1] = P->disp[i] + P->part[i];
      }
      Reorder(P, &T, A); 
      if (!format_out_flag) 
         printf("\nMATRIX REORDERED, time = %e\n",
                 omp_get_wtime() - start);
      if (!format_out_flag){
          printf("\ncolors = %d\n"
                 "partition: mean = %.2f, max %d, min %d\n",
                 P->nparts, MeanInt(P->part, P->nparts),
                 MaxInt(P->part,P->nparts), MinInt(P->part,P->nparts));
          printf("%d & %d & %d & %d & %d\\\\\n\n",
                 nnz, n, P->nparts, MaxInt(P->part,P->nparts),
                 MinInt(P->part,P->nparts));
      }
   }
   else if (block_flag){
      idx_t nparts;
      nparts = (idx_t)P->nparts;
   
      idx_t options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
      idx_t *perm = (idx_t *)calloc(A->n, sizeof(idx_t));
      start = omp_get_wtime();
      flag =  METIS_PartGraphKway(&(G.n), &ncon, G.xadj, G.adjncy, NULL,
                                  NULL, NULL, &nparts, NULL, NULL,
                                  options, &objval, perm);
      stop = omp_get_wtime() - start;
      if (flag != METIS_OK)
         printf("****WARNING****: METIS returned error with code %d.\n", flag);
      if (!format_out_flag)
         printf("\nPartitions determined using METIS, time = %es.\n", stop);
      FreeMetis(&G);
      P->perm = (int *)malloc(A->n * sizeof(int));
      P->part = (int *)calloc(P->nparts, sizeof(int));
      for (int i = 0; i < A->n; i++){
         P->perm[i] = perm[i];
         P->part[P->perm[i]]++;
      }
      P->disp = (int *)malloc((P->nparts+1) * sizeof(int));
      P->disp[0] = 0;
      for (int i = 0; i < P->nparts; i++){
         P->disp[i+1] = P->disp[i] + P->part[i];
      }
      free(perm);
      Reorder(P, &T, A);
      ScaleDiag_CSC(A);
      SEQ_SetupBlocks(*A, *P, D, B);
      *Pard = (PardisoInfo *)malloc(P->nparts * sizeof(PardisoInfo));

      if (loc_direct_flag){
         for (int p = 0; p < P->nparts; p++){
            (*Pard)[p].csr.n = (*D)[p].n;
            (*Pard)[p].csr.nnz = (*D)[p].nnz;
            (*Pard)[p].csr.ja = 
               (MKL_INT *)calloc((*Pard)[p].csr.nnz, sizeof(MKL_INT));
            (*Pard)[p].csr.ia = 
               (MKL_INT *)calloc((*Pard)[p].csr.n+1, sizeof(MKL_INT));
            for (int i = 0; i < (*Pard)[p].csr.nnz; i++){
               (*Pard)[p].csr.ja[i] = (*D)[p].i[i] - P->disp[p];
            }
            for (int i = 0; i < (*Pard)[p].csr.n+1; i++){
               (*Pard)[p].csr.ia[i] = (*D)[p].j_ptr[i];
            }
   
            for (int i = 0; i < 64; i++){
               (*Pard)[p].iparm[i] = 0;
               (*Pard)[p].pt[i] = 0;
            }

          //  for (int i = 0; i < (*Pard)[p].csr.n; i++){
          //     for (int j = 0; j < (*Pard)[p].csr.ia[i+1] - (*Pard)[p].csr.ia[i]; j++){
          //        printf("%d ", (*Pard)[p].csr.ja[(*Pard)[p].csr.ia[i]+j]); 
          //     }
          //     printf("\n");
          //  }
   
            (*Pard)[p].mtype = 11;
            (*Pard)[p].nrhs = 1;
            (*Pard)[p].iparm[17] = -1;
            (*Pard)[p].iparm[18] = -1;
            (*Pard)[p].iparm[0] = 1;         /* No solver default */
            (*Pard)[p].iparm[1] = 0;         /* Fill-in reordering from METIS */
            (*Pard)[p].iparm[7] = 1;
            (*Pard)[p].iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
            (*Pard)[p].iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
            (*Pard)[p].iparm[12] = 1;
            (*Pard)[p].iparm[24] = 1;
            (*Pard)[p].iparm[26] = 1;
            (*Pard)[p].iparm[34] = 1;        /* turn off 1-based indexing */
            (*Pard)[p].maxfct = 1;           /* Maximum number of numerical factorizations. */
            (*Pard)[p].mnum = 1;         /* Which factorization to use. */
            (*Pard)[p].msglvl = 0;           /* Print statistical information in file */
            (*Pard)[p].error = 0;            /* Initialize error flag */
           /* reordering and Symbolic factorization. this step also allocates
            *all memory that is necessary for the factorization. */
            (*Pard)[p].phase = 11;
            PARDISO((*Pard)[p].pt,
                    &((*Pard)[p].maxfct),
                    &((*Pard)[p].mnum),
                    &((*Pard)[p].mtype),
                    &((*Pard)[p].phase),
                    &((*Pard)[p].csr.n),
                    (*D)[p].a,
                    (*Pard)[p].csr.ia,
                    (*Pard)[p].csr.ja,
                   // (*Pard)[p].perm,
                    &((*Pard)[p].idum),
                    &((*Pard)[p].nrhs),
                    (*Pard)[p].iparm,
                    &((*Pard)[p].msglvl),
                    &((*Pard)[p].ddum),
                    &((*Pard)[p].ddum),
                    &((*Pard)[p].error));
           // printf("%d\n", (*Pard)[p].error);
            (*Pard)[p].phase = 22;
            PARDISO((*Pard)[p].pt,
                    &((*Pard)[p].maxfct),
                    &((*Pard)[p].mnum),
                    &((*Pard)[p].mtype),
                    &((*Pard)[p].phase),
                    &((*Pard)[p].csr.n),
                    (*D)[p].a,
                    (*Pard)[p].csr.ia,
                    (*Pard)[p].csr.ja,
                   // (*Pard)[p].perm,
                    &((*Pard)[p].idum),
                    &((*Pard)[p].nrhs),
                    (*Pard)[p].iparm,
                    &((*Pard)[p].msglvl),
                    &((*Pard)[p].ddum),
                    &((*Pard)[p].ddum),
                    &((*Pard)[p].error));
            (*Pard)[p].phase = 33;
         }
      }
   }
  // strcpy(buffer, "metis_matrix_matlab.txt");
  // WriteBlocks_csc(buffer, *D, *B, *P, 1);
   
  // for (int p = 0; p < P->nparts; p++){ 
  //    Write_csc((*D)[p], 1);
  // }

   if (!format_out_flag) printf("\n******************************************************\n\n");
}

void SEQ_SetupNeighb(CSC A, OrderInfo P, SEQ_NeighbInfo *N, int tid)
{
   int i, row, ind, ai;
   int n = A.n;
   int t_low = P.disp[tid];
   int nt = P.part[tid];
   int q;
   N->len = (int *)calloc(nt, sizeof(int));
   N->res = (double **)calloc(nt, sizeof(double *));
   for (int i = 0; i < nt; i++){
      q = 0;
      ai = t_low + i;
      N->len[i] = A.j_ptr[ai+1]-A.j_ptr[ai]-1;
      N->res[i] = (double *)calloc(N->len[i]+1, sizeof(double));
   }
}

void SEQ_SetupBlocks(CSC A,
                     OrderInfo P,
                     CSC **D,
                     CSC **B)
{
   int ind, row, k;
   double elem;

   for (int i = 0; i < A.n; i++){
      QuicksortPair_int_dbl(A.i,
                            A.a,
                            A.j_ptr[i],
                            A.j_ptr[i+1]-1);
   }

   *D = (CSC *)malloc(P.nparts * sizeof(CSC)); 
   *B = (CSC *)malloc(P.nparts * sizeof(CSC));
   
   vector<list<int>> diag_row_list;
   vector<list<int>> offdiag_row_list;
   vector<list<double>> diag_elem_list;
   vector<list<double>> offdiag_elem_list;

   for (int p = 0; p < P.nparts; p++){ 
      diag_row_list.resize(P.part[p]);
      offdiag_row_list.resize(P.part[p]);
      diag_elem_list.resize(P.part[p]);
      offdiag_elem_list.resize(P.part[p]);
      (*D)[p].n = (*B)[p].n = P.part[p];
      for (int i = P.disp[p]; i < P.disp[p+1]; i++){
         for (int j = 0; j < A.j_ptr[i+1] - A.j_ptr[i]; j++){
            ind = A.j_ptr[i]+j;
            row = A.i[ind];
            elem = A.a[ind];
            k = row - P.disp[p];
            if ((k < P.part[p]) && (k >= 0)){
               diag_row_list[i-P.disp[p]].push_back(row);
               diag_elem_list[i-P.disp[p]].push_back(elem);
            }
            else {
               offdiag_row_list[i-P.disp[p]].push_back(row);
               offdiag_elem_list[i-P.disp[p]].push_back(elem);
            }
         }
      }
      List_to_Block(&(*D)[p], diag_row_list, diag_elem_list);
      List_to_Block(&(*B)[p], offdiag_row_list, offdiag_elem_list);
      (*D)[p].diag = (double *)calloc((*D)[p].n, sizeof(double));
      k = 0;
      for (int i = P.disp[p]; i < P.disp[p+1]; i++){
         diag_row_list.clear();
         offdiag_row_list.clear();
         diag_elem_list.clear();
         offdiag_elem_list.clear();
         for (int j = 0; j < A.j_ptr[i+1] - A.j_ptr[i]; j++){
            ind = A.j_ptr[i]+j;
            row = A.i[ind];
            elem = A.a[ind];
            if (i == row){
               (*D)[p].diag[k] = elem;
               k++;
            }
         }
      }
   }
}


void SEQ_SetupNeighbBlocks(OrderInfo P,
                           CSC *B,
                           SEQ_BlockNeighbInfo *Neighb)
{
   int ind, row, k;
   double elem;
   
   
   Neighb->size = (int *)calloc(P.nparts, sizeof(int));
   Neighb->blocks = (int **)calloc(P.nparts, sizeof(int *));
   list<int> neighb_list; 

   for (int p = 0; p < P.nparts; p++){ 
      for (int i = 0; i < B[p].n; i++){
         for (int j = 0; j < B[p].j_ptr[i+1] - B[p].j_ptr[i]; j++){
            ind = B[p].j_ptr[i]+j;
            row = B[p].i[ind];
            for (int pp = 0; pp < P.nparts; pp++){
               if (row < P.disp[pp+1]){
                  if (pp != p){
                     neighb_list.push_back(pp);
                  }
                  break;
               }
            }
         }
      }
      Neighb->size[p] = neighb_list.size();
      Neighb->blocks[p] = (int *)calloc(Neighb->size[p], sizeof(int));
      for (int i = 0; i < Neighb->size[p]; i++){
         Neighb->blocks[p][i] = neighb_list.front();
         neighb_list.pop_front();
      }
      neighb_list.clear(); 
   }
}
