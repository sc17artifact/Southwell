#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Laplace.h"
#include "Multicolor.h"
#include "FileUtils.h"
#include "MatrixUtils.h"
#include "Misc.h"

using namespace std;

void DMEM_ReadMatrix(FILE *in_file,
                     DMEM_ProcInfo *Proc,
                     MetisGraph *G,
                     Triplet *T,
                     int m, int w);
void DMEM_Reorder(OrderInfo *P,
                  DMEM_ProcInfo *Proc,
                  Triplet *T,
                  CSC *A);
void DMEM_DistributeMatrix(CSC A,
                           DMEM_MatrixInfo *Mat,
                           DMEM_ProcInfo *Proc,
                           int disp_glob);
void DMEM_METIS_SetupBoundary(DMEM_MatrixInfo *Mat,
                              DMEM_ProcInfo *Proc);
void DMEM_MC_SetupBoundary(DMEM_MulticolorInfo *MC_info,
                           DMEM_ProcInfo Proc);
void DMEM_METIS_SetupBlocks(DMEM_MatrixInfo *Mat,
                            DMEM_ProcInfo *Proc,
                            vector<list<int>>
                                *unique_bound_ind_list);
void DMEM_MC_SetupBlocks(DMEM_MatrixInfo *Mat,
                         DMEM_ProcInfo *Proc,
                         DMEM_MulticolorInfo *MC_info,
                         int disp_glob,
                         int current_color,
                         int target_color,
                         vector<list<int>>
                            *unique_bound_ind_list);
void DMEM_RMA_METIS_SetupBoundWin_SymMat(DMEM_ProcInfo *Proc,
                                         DMEM_MatrixInfo *Mat,
                                         vector<list<int>>
                                            unique_bound_ind_list);
void DMEM_RMA_CM_SetupBoundWin_SymMat(DMEM_MulticolorInfo *MC_info,
                                      DMEM_ProcInfo Proc,
                                      vector<vector<vector<list<int>>>>
                                         unique_bound_ind_list);
void DMEM_RowPart(OrderInfo *P,
                  int world_size,
                  int extra_rank,
                  int n);

void DMEM_MetisSetup(FILE *in_file,
                     DMEM_MatrixInfo *Mat,
                     DMEM_ProcInfo *Proc,
                     int m,
                     int w)

{
   idx_t ncon = 1, objval, n;
   int flag = METIS_OK;
   int imbalance = 0;
   int point;
   MetisGraph G;
   Triplet T;
   CSC A;

   idx_t nparts;
   nparts = (idx_t)(Proc->world_size + imbalance);
   Mat->P.nparts = Proc->world_size;

   idx_t options[METIS_NOPTIONS];
   char buffer[256];
   FILE *temp_file;
   double start, stop;

   DMEM_ReadMatrix(in_file, Proc, &G, &T, m, w);

   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

   A.n = G.n;
   A.nnz = G.nnz;
   if (Proc->rank == 0){
      idx_t *perm = (idx_t *)calloc(A.n, sizeof(idx_t));
      if (A.n > Proc->world_size){
         start = MPI_Wtime();
         flag =  METIS_PartGraphKway(&(G.n), &ncon, G.xadj, G.adjncy, NULL, 
                                     NULL, NULL, &nparts, NULL, NULL, 
                                     options, &objval, perm);
         stop = MPI_Wtime() - start;
         if (flag != METIS_OK)
            printf("****WARNING****: METIS returned error with code %d.\n", flag);
         if (!format_out_flag) 
            printf("\nPartitions determined using METIS, time = %es.\n", stop);
         FreeMetis(&G);
      }
      else {
         for (int i = 0; i < A.n; i++) perm[i] = i;
      }
      Mat->P.perm = (int *)malloc(A.n * sizeof(int));
      Mat->P.part = (int *)calloc(Mat->P.nparts, sizeof(int));
      for (int i = 0; i < A.n; i++){
         Mat->P.perm[i] = perm[i];
        // printf("%d\n", perm[i]);
         Mat->P.part[Mat->P.perm[i]]++;
      }
      free(perm);
      DMEM_Reorder(&(Mat->P), Proc, &T, &A);
      ScaleDiag_CSC(&A);
      Mat->n_glob = A.n;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   start = MPI_Wtime();
   DMEM_DistributeMatrix(A, Mat, Proc, 0);
   stop = MPI_Wtime() - start;
   if ((Proc->rank == 0) && !format_out_flag)
      printf("\nMatrix distributed, time = %es.\n", stop);
   if (Proc->rank == 0) FreeCSC(&A);
   start = MPI_Wtime();
   DMEM_METIS_SetupBoundary(Mat, Proc);
   stop = MPI_Wtime() - start;
   if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){

      Mat->Pard.csr.n = Mat->n;
      Mat->Pard.csr.nnz = Mat->D.nnz;
      Mat->Pard.csr.ja = (MKL_INT *)calloc(Mat->Pard.csr.nnz, sizeof(MKL_INT));
      Mat->Pard.csr.ia = (MKL_INT *)calloc(Mat->Pard.csr.n+1, sizeof(MKL_INT));
      for (int i = 0; i < Mat->Pard.csr.nnz; i++){
         Mat->Pard.csr.ja[i] = Mat->D.i[i] - Mat->disp;
      }
      for (int i = 0; i < Mat->Pard.csr.n+1; i++){
         Mat->Pard.csr.ia[i] = Mat->D.j_ptr[i];
      }
      
      for (int i = 0; i < 64; i++){
         Mat->Pard.iparm[i] = 0;
         Mat->Pard.pt[i] = 0;
      }
      
      Mat->Pard.mtype = 11;
      Mat->Pard.nrhs = 1;
      Mat->Pard.iparm[17] = -1;
      Mat->Pard.iparm[18] = -1;
      Mat->Pard.iparm[0] = 1;         /* No solver default */
      Mat->Pard.iparm[1] = 0;         /* Fill-in reordering from METIS */
      Mat->Pard.iparm[7] = 1;
      Mat->Pard.iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
      Mat->Pard.iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
      Mat->Pard.iparm[12] = 1;
      Mat->Pard.iparm[24] = 1;
      Mat->Pard.iparm[26] = 1;
      Mat->Pard.iparm[34] = 1;        /* turn off 1-based indexing */
      Mat->Pard.maxfct = 1;           /* Maximum number of numerical factorizations. */
      Mat->Pard.mnum = 1;         /* Which factorization to use. */
      Mat->Pard.msglvl = 0;           /* Print statistical information in file */
      Mat->Pard.error = 0;            /* Initialize error flag */
     /* reordering and Symbolic factorization. this step also allocates
      *all memory that is necessary for the factorization. */
      Mat->Pard.phase = 11;

     // Mat->Pard.perm = (MKL_INT *)calloc(Mat->Pard.csr.n, sizeof(MKL_INT));
     // for (int q = 0; q < Proc->Neighb.size; q++){
     //    for (int i = 0; i < Proc->Bound.num_recv_points[q]; i++){
     //       point = Proc->Bound.recv_points[q][i];
     //       Mat->Pard.perm[point] = 1;
     //    }
     // }
     // Mat->Pard.iparm[30] = 1;
      start = MPI_Wtime();
      PARDISO(Mat->Pard.pt,
              &(Mat->Pard.maxfct),
              &(Mat->Pard.mnum),
              &(Mat->Pard.mtype),
              &(Mat->Pard.phase),
              &(Mat->Pard.csr.n),
              Mat->D.a,
              Mat->Pard.csr.ia,
              Mat->Pard.csr.ja, 
             // Mat->Pard.perm, 
              &(Mat->Pard.idum),
              &(Mat->Pard.nrhs), 
              Mat->Pard.iparm,
              &(Mat->Pard.msglvl),
              &(Mat->Pard.ddum),
              &(Mat->Pard.ddum),
              &(Mat->Pard.error));
      /* numerical factorization */
      Mat->Pard.phase = 22;
      PARDISO(Mat->Pard.pt,
              &(Mat->Pard.maxfct),
              &(Mat->Pard.mnum),
              &(Mat->Pard.mtype),
              &(Mat->Pard.phase),
              &(Mat->Pard.csr.n),
              Mat->D.a,
              Mat->Pard.csr.ia,
              Mat->Pard.csr.ja,
             // Mat->Pard.perm,
              &(Mat->Pard.idum),
              &(Mat->Pard.nrhs),
              Mat->Pard.iparm,
              &(Mat->Pard.msglvl),
              &(Mat->Pard.ddum),
              &(Mat->Pard.ddum),
              &(Mat->Pard.error));
      Mat->Pard.wtime_setup = MPI_Wtime() - start;
      Mat->Pard.phase = 33;
   }
   FreeCSC(&(Mat->A));
   if ((Proc->rank == 0) && !format_out_flag)
      printf("\nBoundary point information initialized, time = %es.\n", stop);
  // strcpy(buffer, "metis_matrix_matlab.txt");
  // if (Proc->rank == 0) remove(buffer);
  // MPI_Barrier(MPI_COMM_WORLD);
  // DMEM_WriteBlocks_csc(buffer, *Mat, *Proc, 1);
  // for (int p = 0; p < Proc->world_size; p++){
  //    if (Proc->rank == p){
  //       printf("rank %d, n = %d:\n\t", p, Mat->n);
  //       for (int q = 0; q < Proc->Neighb.size; q++){
  //          printf("%d, %d:\n\t\t", Proc->Neighb.ranks[q], Proc->Bound.Rma.recv_count);
  //          for (int i = 0; i < Proc->Bound.num_send_points[q]; i++){
  //             printf("%d ", Proc->Bound.send_points[q][i]);
  //          }
  //          printf("\n\t");
  //       }
  //       printf("\n");
  //    } 
  //    MPI_Barrier(MPI_COMM_WORLD);
  // }
  // Write_csc(Mat->A, 1);
   if ((Proc->rank == 0) && !format_out_flag) 
      printf("\n|******************************************|\n\n");
}

void DMEM_MulticolorSetup(FILE *in_file,
                          DMEM_MulticolorInfo *MC_info,
                          DMEM_ProcInfo *Proc,
                          int m,
                          int w)

{
   idx_t ncon = 1, objval, n;
   int flag = METIS_OK;
   int imbalance = 0;
   MetisGraph G;
   Triplet T;
   CSC A;

   char buffer[256];
   FILE *temp_file;
   double start, stop;

   DMEM_ReadMatrix(in_file, Proc, &G, &T, m, w);
   if (Proc->rank == 0){
      Metis_to_CSC(&A, G);
      FreeMetis(&G);
      MC_info->P_glob.perm = (int *)malloc(A.n * sizeof(int));
      if (!format_out_flag) printf("\nREORDERING USING MULTICOLORING\n");
      Multicolor(A, &(MC_info->P_glob));
      if (!format_out_flag){
         printf("colors = %d\n"
                "partition: mean = %.2f, max %d, min %d\n",
                MC_info->P_glob.nparts, MeanInt(MC_info->P_glob.part, MC_info->P_glob.nparts),
                MaxInt(MC_info->P_glob.part, MC_info->P_glob.nparts),
                MinInt(MC_info->P_glob.part, MC_info->P_glob.nparts));
         printf("%d & %d & %d & %d & %d\\\\\n",
                A.nnz, A.n, MC_info->P_glob.nparts,
                MaxInt(MC_info->P_glob.part, MC_info->P_glob.nparts),
                MinInt(MC_info->P_glob.part, MC_info->P_glob.nparts));
      }
      DMEM_Reorder(&(MC_info->P_glob), Proc, &T, &A);
      MC_info->n_glob = A.n;
      MC_info->nnz_glob = A.nnz;
     // Write_csc(A,1);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&(MC_info->n_glob), 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&(MC_info->nnz_glob), 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&(MC_info->P_glob.nparts), 1, MPI_INT, 0, MPI_COMM_WORLD); 
   if (Proc->rank > 0){
      MC_info->P_glob.part = (int *)calloc(MC_info->P_glob.nparts, sizeof(int));
      MC_info->P_glob.disp = (int *)calloc(MC_info->P_glob.nparts+1, sizeof(int));
   }
   MC_info->Mat = (DMEM_MatrixInfo *)malloc(MC_info->P_glob.nparts * sizeof(DMEM_MatrixInfo));
   MC_info->B = (CSC ***)malloc(MC_info->P_glob.nparts * sizeof(CSC **));
   for (int c = 0; c < MC_info->P_glob.nparts; c++)
      MC_info->B[c] = (CSC **)malloc(MC_info->P_glob.nparts * sizeof(CSC *));
   MC_info->NeighbSend = 
      (DMEM_NeighbInfo **)malloc(MC_info->P_glob.nparts * sizeof(DMEM_NeighbInfo *));
   for (int c = 0; c < MC_info->P_glob.nparts; c++)
      MC_info->NeighbSend[c] =
         (DMEM_NeighbInfo *)malloc(MC_info->P_glob.nparts * sizeof(DMEM_NeighbInfo));
   MC_info->NeighbRecv =
      (DMEM_NeighbInfo **)malloc(MC_info->P_glob.nparts * sizeof(DMEM_NeighbInfo *));
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      MC_info->NeighbRecv[c] =
         (DMEM_NeighbInfo *)malloc(MC_info->P_glob.nparts * sizeof(DMEM_NeighbInfo));
   }
   MC_info->Bound =
      (DMEM_BoundInfo  **)malloc(MC_info->P_glob.nparts * sizeof(DMEM_BoundInfo *));
   for (int c = 0; c < MC_info->P_glob.nparts; c++)
      MC_info->Bound[c] =
         (DMEM_BoundInfo *)malloc(MC_info->P_glob.nparts * sizeof(DMEM_BoundInfo));

   MC_info->D = (CSC *)malloc(MC_info->P_glob.nparts * sizeof(CSC));
   MC_info->disp_loc = (int *)calloc((MC_info->P_glob.nparts+1), sizeof(int));
   MC_info->P_loc = (OrderInfo *)malloc(MC_info->P_glob.nparts * sizeof(OrderInfo));
   MPI_Barrier(MPI_COMM_WORLD);
   MC_info->n = 0;
   if ((Proc->rank == 0) && !format_out_flag)
         printf("\n");
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      MC_info->B[c] = (CSC **)malloc(MC_info->P_glob.nparts * sizeof(CSC *));
      MPI_Bcast(&(MC_info->P_glob.disp[c]), 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&(MC_info->P_glob.part[c]), 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (Proc->rank == 0){
         DMEM_RowPart(&(MC_info->P_loc[c]), 
                      Proc->world_size,
                      Proc->world_size-1, 
                      MC_info->P_glob.part[c]);
         MC_info->Mat[c].P.part = (int *)calloc(Proc->world_size, sizeof(int));
         MC_info->Mat[c].P.disp = (int *)calloc(Proc->world_size+1, sizeof(int));
         for (int i = 0; i < Proc->world_size; i++){
            MC_info->Mat[c].P.part[i] = MC_info->P_loc[c].part[i];   
            MC_info->Mat[c].P.disp[i+1] = MC_info->P_loc[c].disp[i+1];
         }
         MC_info->Mat[c].n_glob = MC_info->P_glob.part[c];
      }
      else {
         MC_info->P_loc[c].part = (int *)calloc(Proc->world_size, sizeof(int));
         MC_info->P_loc[c].disp = (int *)calloc(Proc->world_size+1, sizeof(int));
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(MC_info->P_loc[c].part, 
                Proc->world_size, 
                MPI_INT, 
                0, 
                MPI_COMM_WORLD);
      MPI_Bcast(MC_info->P_loc[c].disp, 
                Proc->world_size+1, 
                MPI_INT, 
                0, 
                MPI_COMM_WORLD);
      MC_info->P_loc[c].nparts = Proc->world_size;
      MPI_Barrier(MPI_COMM_WORLD);
      MC_info->Mat[c].P.nparts = Proc->world_size;
      start = MPI_Wtime();
      DMEM_DistributeMatrix(A, 
                            &(MC_info->Mat[c]), 
                            Proc, 
                            MC_info->P_glob.disp[c]);
      stop = MPI_Wtime() - start;
      MC_info->D[c].n = MC_info->Mat[c].A.n;
      MC_info->n += MC_info->Mat[c].P.part[Proc->rank];
      MC_info->disp_loc[c+1] = 
         MC_info->disp_loc[c] + MC_info->P_loc[c].part[Proc->rank];
      MPI_Barrier(MPI_COMM_WORLD);
      if ((Proc->rank == 0) && !format_out_flag)
         printf("Matrix distributed for color %d, time = %e\n",
                c, MPI_Wtime()-start);
   }
   MC_info->P_glob.disp[MC_info->P_glob.nparts] = MC_info->n_glob;
   if (Proc->rank == 0) FreeCSC(&A);
   MPI_Barrier(MPI_COMM_WORLD); 
   start = MPI_Wtime();
   DMEM_MC_SetupBoundary(MC_info, *Proc);
   MPI_Barrier(MPI_COMM_WORLD);
   if ((Proc->rank == 0) && !format_out_flag)
      printf("Neighbor information intialized, time = %e\n",
              MPI_Wtime()-start);
  // strcpy(buffer, "mc_matrix_matlab.txt");
  // if (Proc->rank == 0) remove(buffer);
  // MPI_Barrier(MPI_COMM_WORLD);
  // DMEM_WriteBlocks_MC(buffer, 
  //                     *MC_info, 
  //                     *Proc, 
  //                     1);
   if ((Proc->rank == 0) && !format_out_flag)
      printf("\n|******************************************|\n\n");
   free(MC_info->Mat);
}

void DMEM_ReadMatrix(FILE *in_file,
                     DMEM_ProcInfo *Proc,
                     MetisGraph *G, 
                     Triplet *T, 
                     int m, int w)
{
   double start, stop;
   if ((Proc->rank == 0) && (!format_out_flag))
      printf("\n|********** INITIALIZING PROBLEM **********|\n");
   start = MPI_Wtime();
   if (mat_file_flag){
      if (Proc->rank == 0){
         if (!format_out_flag){
            printf("\nReading matrix from file.\n");
         }
         ReadBinary_fread_metis(in_file, G, T);
      }
   }
   else {
      if (Proc->rank == 0){
         if (!format_out_flag){
            printf("\nUsing 5-pt centered difference discretization of the Laplace PDE.\n");
         }
         Laplace2D_FD5pt_metis(G, T, m, w);
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   stop = MPI_Wtime() - start;
   if ((Proc->rank == 0) && !format_out_flag){
      printf("\nMatrix loaded, time = %es, nnz = %d, n = %d.\n",
             stop, G->nnz, G->n);
   }
}

void DMEM_Reorder(OrderInfo *P, 
                  DMEM_ProcInfo *Proc,
                  Triplet *T,
                  CSC *A)
{
   double start, stop;
   P->disp = (int *)malloc((P->nparts+1) * sizeof(int));
   P->disp[0] = 0;
   for (int i = 0; i < P->nparts; i++){
      P->disp[i+1] = P->disp[i] + P->part[i];
   }
   start = MPI_Wtime();
   ReorderTriplet_csc(*T, A, P);
   stop = MPI_Wtime() - start;
   free(T->a);
   free(T->i);
   free(T->j);
   free(P->perm);
   free(P->map);
   if (!format_out_flag)
      printf("\nMatrix reordered, time = %es.\n", stop);
         // Write_csc(A, 1);
}

void DMEM_DistributeMatrix(CSC A,
                           DMEM_MatrixInfo *Mat, 
                           DMEM_ProcInfo *Proc,
                           int disp_glob)
{
   int disp;
   int j_ptr_disp;
   Mat->P.dispv = (int *)malloc(Mat->P.nparts * sizeof(int));
   int *j_ptr_extra = (int *)malloc(Mat->P.nparts * sizeof(int));
   if (Proc->rank == 0){
      for (int i = 0; i < Mat->P.nparts; i++){ 
         Mat->P.dispv[i] = Mat->P.disp[i];
         j_ptr_extra[i] = A.j_ptr[disp_glob+Mat->P.disp[i+1]];
      }
   }
   else {
      Mat->P.part = (int *)malloc(Mat->P.nparts * sizeof(int));
      Mat->P.disp = (int *)malloc((Mat->P.nparts+1) * sizeof(int));
   }
   
   MPI_Bcast(Mat->P.part, Mat->P.nparts, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(Mat->P.disp, Mat->P.nparts+1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(Mat->P.dispv, Mat->P.nparts, MPI_INT, 0, MPI_COMM_WORLD);
   
   MPI_Bcast(&(Mat->n_glob), 1, MPI_INT, 0, MPI_COMM_WORLD);
  // MPI_Bcast(&(Mat->nnz), 1, MPI_INT, 0, MPI_COMM_WORLD);   

   Mat->disp = Mat->P.disp[Proc->rank];
   Mat->n = Mat->A.n = Mat->D.n = Mat->P.part[Proc->rank];
   Mat->A.diag = (double *)calloc(Mat->n, sizeof(double)); 
   MPI_Scatterv(&(A.diag[disp_glob]), Mat->P.part, Mat->P.dispv, MPI_DOUBLE,
                Mat->A.diag, Mat->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   Mat->A.j_ptr = (int *)calloc(Mat->n+1, sizeof(int));
   MPI_Scatterv(&(A.j_ptr[disp_glob]), 
                Mat->P.part, 
                Mat->P.dispv, 
                MPI_INT, 
                Mat->A.j_ptr, 
                Mat->n+1, 
                MPI_INT, 
                0, 
                MPI_COMM_WORLD);
   MPI_Scatter(j_ptr_extra, 
               1, 
               MPI_INT, 
               &(Mat->A.j_ptr[Mat->n]), 
               1, 
               MPI_INT,
               0, 
               MPI_COMM_WORLD);
   j_ptr_disp = Mat->A.j_ptr[0];
   for (int i = 0; i < (Mat->n+1); i++){
      Mat->A.j_ptr[i] -= j_ptr_disp; 
   }
   
   int *dispv_nnz = (int *)malloc(Mat->P.nparts * sizeof(int)); 
   int *part_nnz = (int *)malloc(Mat->P.nparts * sizeof(int));
   int count_nnz, count_nnz_prev;
   int tot_nnz = 0;
   if (Proc->rank == 0){
      dispv_nnz[0] = 0;
      for (int i = 0; i < Mat->P.nparts; i++){
         count_nnz = 0;
         disp = Mat->P.disp[i];
         for (int j = disp; j < (disp + Mat->P.part[i]); j++){
            count_nnz += (A.j_ptr[disp_glob+j+1] - A.j_ptr[disp_glob+j]);
         }
         if (i > 0) dispv_nnz[i] = dispv_nnz[i-1] + count_nnz_prev;
         part_nnz[i] = count_nnz_prev = count_nnz;
      }
   }
   
  // printf("%d, %d\n", A.n, disp_glob);

   MPI_Bcast(part_nnz, Mat->P.nparts, MPI_INT, 0, MPI_COMM_WORLD);
   Mat->A.nnz = Mat->nnz = part_nnz[Proc->rank];
   Mat->A.i = (int *)calloc(Mat->A.nnz, sizeof(int));
   Mat->A.a = (double *)calloc(Mat->A.nnz, sizeof(double));

   if (Proc->rank == 0){
      disp = A.j_ptr[disp_glob];
   }
   MPI_Bcast(&disp, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&(A.i[disp]), part_nnz, dispv_nnz, MPI_INT,
                Mat->A.i, Mat->A.nnz, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(&(A.a[disp]), part_nnz, dispv_nnz, MPI_DOUBLE,
                Mat->A.a, Mat->A.nnz, MPI_DOUBLE, 0,MPI_COMM_WORLD);

   free(part_nnz);
   free(dispv_nnz);
   free(j_ptr_extra);
}

void DMEM_METIS_SetupBoundary(DMEM_MatrixInfo *Mat,
                              DMEM_ProcInfo *Proc)
{
   int cum_part, next_part;
   int row_rank;
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   int break_color_loop;
   double elem;
   vector<list<int>> unique_bound_ind_list(Proc->world_size);
   DMEM_METIS_SetupBlocks(Mat, Proc, &(unique_bound_ind_list));

   Proc->Res.norm =
      (double *)calloc(Proc->Neighb.size+1, sizeof(double));
   Proc->Bound.num_send_points =
       (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.row_min =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.row_max =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.send_points =
      (int **)calloc(Proc->Neighb.size, sizeof(int *));
   Proc->Bound.num_send_points =
      (int *)calloc(Proc->world_size, sizeof(int));
   Proc->Bound.row_range_minmax =
      (int *)calloc(Proc->world_size, sizeof(int));
   Proc->Bound.row_range_disp =
      (int *)calloc(Proc->world_size, sizeof(int));
   Proc->Bound.Rma.send_count =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.send_disp =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.num_recv_points =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.recv_disp =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.send =
      (double **)calloc(Proc->Neighb.size, sizeof(double *));
   Proc->Bound.recv_points =
      (int **)calloc(Proc->Neighb.size, sizeof(int *));
   Proc->Bound.points =
         (int **)calloc(Proc->Neighb.size, sizeof(int *));

   if (solver_flag == DMEM_SOLVER_DS){
      Proc->Res.norm =
         (double *)calloc(Proc->Neighb.size+1, sizeof(double));
      Proc->Res.norm_prev =
         (double *)calloc(Proc->Neighb.size+1, sizeof(double));
      Proc->Res.neighb_explicit_flag =
         (int *)calloc(Proc->Neighb.size, sizeof(int));
      Proc->Bound.Rma.send_prev =
         (double **)calloc(Proc->Neighb.size, sizeof(double *));
      Proc->Res.neighb_norm_estim =
         (double *)calloc(Proc->Neighb.size, sizeof(double));
      Proc->Res.norm_estim_squ =
         (double *)calloc(Proc->Neighb.size, sizeof(double));
      Proc->Res.my_norm_estim_squ =
         (double *)calloc(Proc->Neighb.size, sizeof(double));
      Proc->Res.my_norm_estim_squ_prev =
         (double *)calloc(Proc->Neighb.size, sizeof(double));
      Proc->Bound.res =
         (double **)calloc(Proc->Neighb.size, sizeof(double *));

   }
   else if (solver_flag == DMEM_SOLVER_PS){
      Proc->Res.norm_prev =
         (double *)calloc(Proc->Neighb.size+1, sizeof(double));
      Proc->Res.norm_send_prev =
         (double *)calloc(Proc->Neighb.size+1, sizeof(double));
   }
  
   DMEM_RMA_METIS_SetupBoundWin_SymMat(Proc, 
                                       Mat, 
                                       unique_bound_ind_list);
}

void DMEM_METIS_SetupBlocks(DMEM_MatrixInfo *Mat,
                            DMEM_ProcInfo *Proc,
                            vector<list<int>>
                                *unique_bound_ind_list)
{
  
   int cum_part, next_part;
   int row_rank;
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   int break_color_loop;
   double elem; 

   /* Sort columns of un-split matrix (required by PARDISO). */
   for (int i = 0; i < Mat->n; i++){
      QuicksortPair_int_dbl(Mat->A.i, 
                            Mat->A.a, 
                            Mat->A.j_ptr[i], 
                            Mat->A.j_ptr[i+1]-1);
   }

   
   list<int> neighb_list;
   vector<list<int>> diag_ind_list(Mat->n);
   vector<list<double>> diag_elem_list(Mat->n);
   vector<list<int>> bound_ind_list(Proc->world_size);
   vector<vector<list<int>>> 
      bound_ind_list_2d(Proc->world_size, 
                        vector<list<int>>(Mat->n));
   vector<vector<list<double>>>
      bound_elem_list_2d(Proc->world_size,
                         vector<list<double>>(Mat->n));
   for (int i = 0; i < Mat->n; i++){
      for (int j = 0; j < Mat->A.j_ptr[i+1]-Mat->A.j_ptr[i]; j++){
         k = Mat->A.j_ptr[i]+j;
         row = Mat->A.i[k];
         elem = Mat->A.a[k];
         row_rank = 0;
         while (Mat->P.disp[row_rank] < Mat->n_glob){
            next_part = Mat->P.part[row_rank];
            if (row < (Mat->P.disp[row_rank] + next_part)){
               if (row_rank == Proc->rank){
                  diag_ind_list[i].push_back(row);
                  diag_elem_list[i].push_back(elem);
               }
               else {
                  neighb_list.push_back(row_rank);
                  neighb_list.sort();
                  neighb_list.unique();
                  bound_ind_list[row_rank].push_back(row);

                  bound_ind_list_2d[row_rank][i].push_back(row);
                  bound_elem_list_2d[row_rank][i].push_back(elem);
               }
               break;
            }
            row_rank++;
         }
      }
   }
   Proc->Neighb.size = neighb_list.size();
   Proc->Neighb.ranks = (int *)calloc(Proc->Neighb.size, sizeof(int *));
   for (int i = 0; i < Proc->Neighb.size; i++){
      neighb_rank = neighb_list.front();

      neighb_list.pop_front();
      Proc->Neighb.ranks[i] = neighb_rank;

      bound_ind_list[neighb_rank].sort();
      (*unique_bound_ind_list)[neighb_rank].assign(
         bound_ind_list[neighb_rank].begin(), 
         bound_ind_list[neighb_rank].end());
      (*unique_bound_ind_list)[neighb_rank].unique();
   }

   List_to_Block(&(Mat->D), diag_ind_list, diag_elem_list);

   Mat->D.diag = (double *)calloc(Mat->n, sizeof(double));
   k = 0;
   for (int i = 0; i < Mat->n; i++){
      for (int j = 0; j < Mat->D.j_ptr[i+1]-Mat->D.j_ptr[i]; j++){
         ind = Mat->D.j_ptr[i]+j;
         row = Mat->D.i[ind];
         elem = Mat->D.a[ind];
         if ((Mat->disp+i) == row){
            Mat->D.diag[k] = elem;
            k++;
         }
      }
   }
  

   Mat->B = (CSC *)malloc(Proc->Neighb.size * sizeof(CSC));
   for (int i = 0; i < Proc->Neighb.size; i++){
      Mat->B[i].n = Mat->n;
      neighb_rank = Proc->Neighb.ranks[i];
      List_to_Block(&(Mat->B[i]), bound_ind_list_2d[neighb_rank], 
                    bound_elem_list_2d[neighb_rank]);
   }
}


void DMEM_RMA_METIS_SetupBoundWin_SymMat(DMEM_ProcInfo *Proc,
                                         DMEM_MatrixInfo *Mat,
                                         vector<list<int>>
                                            unique_bound_ind_list)
{
   int row, ind, shift_row;
   int neighb_rank;
   int k, h;
   int map;
   int *sendbuff, *recvbuff;
   int *sendcnts, *recvcnts;
   int *sdispls, *rdispls;
   double elem;
   MPI_Request *req_send, req_recv;
   MPI_Status *stat_send, stat_recv;

   req_send = (MPI_Request *)malloc(Proc->Neighb.size * sizeof(MPI_Request));
   stat_send = (MPI_Status *)malloc(Proc->Neighb.size * sizeof(MPI_Status));

   Proc->Bound.Rma.recv_count = 0;
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      Proc->Bound.num_send_points[q] = 
         unique_bound_ind_list[neighb_rank].size();
      Proc->Bound.row_min[q] = unique_bound_ind_list[neighb_rank].front();
      Proc->Bound.row_max[q] = unique_bound_ind_list[neighb_rank].back();
      Proc->Bound.row_range_minmax[q] =
         (Proc->Bound.row_max[q] - Proc->Bound.row_min[q] + 1);
      Proc->Bound.row_range_disp[q] =
         (Proc->Bound.row_max[q] - Mat->P.disp[neighb_rank] + 1);
      Proc->Bound.Rma.send_count[q] = Proc->Bound.num_send_points[q];
      if (solver_flag == DMEM_SOLVER_PS){
         Proc->Bound.Rma.send_count[q]++;
      }
      Proc->Bound.points[q] = 
         (int *)calloc(Proc->Bound.num_send_points[q], sizeof(int));
      for (int i = 0; i < Proc->Bound.num_send_points[q]; i++){
         Proc->Bound.points[q][i] = unique_bound_ind_list[neighb_rank].front();
         unique_bound_ind_list[neighb_rank].pop_front();
      }
      
      Proc->Bound.send_points[q] =
         (int *)calloc(Proc->Bound.row_range_minmax[q], sizeof(int));
      for (int i = 0; i < Proc->Bound.row_range_minmax[q]; i++)
         Proc->Bound.send_points[q][i] = -1;
      for (int i = 0; i < Proc->Bound.num_send_points[q]; i++){
         ind = Proc->Bound.points[q][i] - Proc->Bound.row_min[q];
         Proc->Bound.send_points[q][ind] = i;
         Proc->Bound.points[q][i] -= Mat->P.disp[neighb_rank];
      }
   }

   sendbuff = (int *)calloc(Proc->world_size, sizeof(int));
   recvbuff = (int *)calloc(Proc->world_size, sizeof(int));
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      sendbuff[neighb_rank] = Proc->Bound.num_send_points[q];
   }
   MPI_Alltoall(sendbuff, 1, MPI_INT, recvbuff, 1, MPI_INT, MPI_COMM_WORLD); 
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      Proc->Bound.num_recv_points[q] = recvbuff[neighb_rank];
      if (solver_flag == DMEM_SOLVER_DS){
         Proc->Bound.Rma.send_count[q] += (Proc->Bound.num_recv_points[q] + 2);
         Proc->Bound.Rma.send_prev[q] =
            (double *)calloc(Proc->Bound.Rma.send_count[q], sizeof(double));
      }
      Proc->Bound.Rma.send[q] =
         (double *)calloc(Proc->Bound.Rma.send_count[q], sizeof(double));

      Proc->Bound.Rma.recv_disp[q] = Proc->Bound.Rma.recv_count;
      Proc->Bound.Rma.recv_count += Proc->Bound.num_recv_points[q];
      if (solver_flag == DMEM_SOLVER_DS){
         Proc->Bound.Rma.recv_count +=
            (Proc->Bound.num_send_points[q] + 2);
         Proc->Bound.res[q] =
            (double *)calloc(Proc->Bound.num_send_points[q],
                                            sizeof(double));
      }
      else if (solver_flag == DMEM_SOLVER_PS){
         Proc->Bound.Rma.recv_count++;
      }
      Proc->Bound.recv_points[q] =
         (int *)calloc(Proc->Bound.num_recv_points[q], sizeof(int));
   }
   free(sendbuff);
   free(recvbuff);

   sendbuff = (int *)calloc(Proc->world_size, sizeof(int));
   recvbuff = (int *)calloc(Proc->world_size, sizeof(int));
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      sendbuff[neighb_rank] = Proc->Bound.Rma.recv_disp[q];
   }
   MPI_Alltoall(sendbuff, 1, MPI_INT, recvbuff, 1, MPI_INT, MPI_COMM_WORLD);
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      Proc->Bound.Rma.send_disp[q] = recvbuff[neighb_rank];
   }
   free(sendbuff);
   free(recvbuff);

   sendcnts = (int *)calloc(Proc->world_size, sizeof(int));
   recvcnts = (int *)calloc(Proc->world_size, sizeof(int));
   sdispls  = (int *)calloc((Proc->world_size+1), sizeof(int));
   rdispls  = (int *)calloc((Proc->world_size+1), sizeof(int));

   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      sendcnts[neighb_rank] = Proc->Bound.num_send_points[q];
      recvcnts[neighb_rank] = Proc->Bound.num_recv_points[q];
   }
   for (int p = 0; p < Proc->world_size; p++){
      sdispls[p+1] = sdispls[p] + sendcnts[p];
      rdispls[p+1] = rdispls[p] + recvcnts[p];
   }

   sendbuff = (int *)calloc(sdispls[Proc->world_size], sizeof(int));
   recvbuff = (int *)calloc(rdispls[Proc->world_size], sizeof(int));

   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      for (int i = 0; i < Proc->Bound.num_send_points[q]; i++){
         sendbuff[sdispls[neighb_rank]+i] 
            = Proc->Bound.points[q][i];
      }
   }
   MPI_Alltoallv(sendbuff,
                 sendcnts,
                 sdispls,
                 MPI_INT,
                 recvbuff,
                 recvcnts,
                 rdispls,
                 MPI_INT,
                 MPI_COMM_WORLD);
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      for (int i = 0; i < Proc->Bound.num_recv_points[q]; i++){
         Proc->Bound.recv_points[q][i] = recvbuff[rdispls[neighb_rank]+i];
      }
   }

   free(sendbuff);
   free(recvbuff);
   free(sendcnts);
   free(recvcnts);
   free(sdispls);
   free(rdispls);

   MPI_Win_allocate(sizeof(double) * Proc->Bound.Rma.recv_count,
                    sizeof(double),
                    MPI_INFO_NULL,
                    MPI_COMM_WORLD,
                    &(Proc->Bound.Rma.recv),
                    &(Proc->Bound.Rma.win));
   Proc->Bound.Rma.recv_prev =
      (double *)calloc(Proc->Bound.Rma.recv_count, sizeof(double));
   for (int i = 0; i < Proc->Bound.Rma.recv_count; i++){
      Proc->Bound.Rma.recv[i] = 0;
   }
   if ((solver_flag == DMEM_SOLVER_PS) && (implement_flag == DMEM_SOS)){
         Proc->Res.Rma.send_disp =
            (int *)calloc(Proc->Bound.Rma.recv_count, sizeof(int));
         sendbuff = (int *)calloc(Proc->world_size, sizeof(int));
         recvbuff = (int *)calloc(Proc->world_size, sizeof(int));
         for (int q = 0; q < Proc->Neighb.size; q++){
            neighb_rank = Proc->Neighb.ranks[q];
            sendbuff[neighb_rank] = q;
         }
         MPI_Alltoall(sendbuff, 1, MPI_INT, recvbuff, 1, MPI_INT, MPI_COMM_WORLD);
         for (int q = 0; q < Proc->Neighb.size; q++){
            neighb_rank = Proc->Neighb.ranks[q];
            Proc->Res.Rma.send_disp[q] = recvbuff[neighb_rank];
         }
         free(sendbuff);
         free(recvbuff);
         MPI_Win_allocate(sizeof(double) * Proc->Neighb.size,
                          sizeof(double),
                          MPI_INFO_NULL,
                          MPI_COMM_WORLD,
                          &(Proc->Res.Rma.recv),
                          &(Proc->Res.Rma.win));
         Proc->Res.Rma.recv_prev =
            (double *)calloc(Proc->Neighb.size, sizeof(double));
         for (int q = 0; q < Proc->Neighb.size; q++){
            Proc->Res.Rma.recv[q] = 0;
         }
   }

   Proc->Conv.win_size = 3;
   MPI_Win_allocate(sizeof(int) * Proc->Conv.win_size,
                    sizeof(int),
                    MPI_INFO_NULL,
                    MPI_COMM_WORLD,
                    &(Proc->Conv.recv),
                    &(Proc->Conv.win));
   for (int i = 0; i < Proc->Conv.win_size; i++){
      Proc->Conv.recv[i] = 0;
   }
   Proc->Conv.my_flag = 0;
   Proc->Conv.phase = 1;

   MPI_Group world_group;
   MPI_Comm_group(MPI_COMM_WORLD, &world_group);
   MPI_Group_incl(world_group, Proc->Neighb.size,
                  Proc->Neighb.ranks, &(Proc->Neighb.mpi_group));

}

void DMEM_MC_SetupBoundary(DMEM_MulticolorInfo *MC_info,
                           DMEM_ProcInfo Proc)
{
   int cum_part, next_part;
   int row_rank;
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   int break_color_loop;
   int *temp_send, *temp_recv, temp_size;
   double elem;

   vector<vector<vector<list<int>>>>
      unique_bound_ind_list(MC_info->P_glob.nparts,
         vector<vector<list<int>>>(MC_info->P_glob.nparts,
            vector<list<int>>(Proc.world_size)));

   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
         DMEM_MC_SetupBlocks(&(MC_info->Mat[c]),
                             &Proc,
                             MC_info,
                             MC_info->P_glob.disp[c],
                             c,
                             cc,
                             &(unique_bound_ind_list[c][cc]));
      }
      FreeCSC(&(MC_info->Mat[c].A));
   }

   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
         MC_info->Bound[c][cc].num_send_points =
             (int *)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].row_min =
            (int *)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].row_max =
            (int *)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].send_points =
            (int **)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int *));
         MC_info->Bound[c][cc].num_send_points =
            (int *)calloc(Proc.world_size, sizeof(int));
         MC_info->Bound[c][cc].row_range_minmax =
            (int *)calloc(Proc.world_size, sizeof(int));
         MC_info->Bound[c][cc].row_range_disp =
            (int *)calloc(Proc.world_size, sizeof(int));
         MC_info->Bound[c][cc].Rma.send_count =
            (int *)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].Rma.send_disp =
            (int *)calloc(MC_info->NeighbSend[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].Rma.send =
            (double **)calloc(MC_info->NeighbSend[c][cc].size, 
                                            sizeof(double *));
         MC_info->Bound[c][cc].points =
               (int **)calloc(MC_info->NeighbSend[c][cc].size, 
                                               sizeof(int *));
      }
   }

   DMEM_RMA_CM_SetupBoundWin_SymMat(MC_info,
                                    Proc,
                                    unique_bound_ind_list);
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            free(MC_info->Bound[c][cc].points[q]);
         }
         free(MC_info->Bound[c][cc].points);
      }
   }
}

void DMEM_MC_SetupBlocks(DMEM_MatrixInfo *Mat,
                         DMEM_ProcInfo *Proc,
                         DMEM_MulticolorInfo *MC_info,
                         int disp_glob,
                         int current_color,
                         int target_color,
                         vector<list<int>>
                            *unique_bound_ind_list)
{
  
   int cum_part, next_part;
   int row_rank;
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   int break_color_loop;
   int diag_flag = 0, offd_flag = 0;
   double elem; 

   /* Sort columns of un-split matrix (required by PARDISO). */
   for (int i = 0; i < Mat->n; i++){
      QuicksortPair_int_dbl(Mat->A.i, 
                            Mat->A.a, 
                            Mat->A.j_ptr[i], 
                            Mat->A.j_ptr[i+1]-1);
   }

   
   list<int> neighb_list;
   vector<list<int>> diag_ind_list(Mat->n);
   vector<list<double>> diag_elem_list(Mat->n);
   vector<list<int>> bound_ind_list(Proc->world_size);
   vector<vector<list<int>>> 
      bound_ind_list_2d(Proc->world_size, 
                        vector<list<int>>(Mat->n));
   vector<vector<list<double>>>
      bound_elem_list_2d(Proc->world_size,
                         vector<list<double>>(Mat->n));
   for (int i = 0; i < Mat->n; i++){
      for (int j = 0; j < Mat->A.j_ptr[i+1]-Mat->A.j_ptr[i]; j++){
         k = Mat->A.j_ptr[i]+j;
         row = Mat->A.i[k];
         elem = Mat->A.a[k];
         break_color_loop = 0;
         for (int c = 0; c < MC_info->P_glob.nparts; c++){
            row_rank = 0;
            while (MC_info->P_loc[c].disp[row_rank] < 
                   MC_info->Mat[c].n_glob){
               next_part = MC_info->Mat[c].P.part[row_rank];
               if (row < (MC_info->P_glob.disp[c] + 
                          MC_info->Mat[c].P.disp[row_rank] + 
                          next_part)){
                  if ((row_rank == Proc->rank) && 
                      (c == current_color)){
                     diag_ind_list[i].push_back(row);
                     diag_elem_list[i].push_back(elem);
                     diag_flag = 1;
                  }
                  else if (c == target_color){
                     neighb_list.push_back(row_rank);
                     neighb_list.sort();
                     neighb_list.unique();
                     bound_ind_list[row_rank].push_back(row);

                     bound_ind_list_2d[row_rank][i].push_back(row);
                     bound_elem_list_2d[row_rank][i].push_back(elem);
                     offd_flag = 1;
                  }
                  break_color_loop = 1;
                  break;
               }
               row_rank++;
            }
            if (break_color_loop) break;
         }
         
      }
   }


   if (diag_flag){
      List_to_Block(&(MC_info->D[current_color]), 
                    diag_ind_list, 
                    diag_elem_list); 
      MC_info->D[current_color].diag = 
         (double *)calloc(MC_info->D[current_color].n, sizeof(double));
      k = 0;
      for (int i = 0; i < MC_info->D[current_color].n; i++){
         for (int j = 0; j < 
                 MC_info->D[current_color].j_ptr[i+1] - 
                    MC_info->D[current_color].j_ptr[i]; j++){
            ind = MC_info->D[current_color].j_ptr[i]+j;
            row = MC_info->D[current_color].i[ind];
            elem = MC_info->D[current_color].a[ind];
            if ((MC_info->P_glob.disp[current_color] +
                 MC_info->P_loc[current_color].disp[Proc->rank]+i) == row){
                MC_info->D[current_color].diag[k] = elem;
               k++;
            }
         }
      }
   }
   MC_info->NeighbSend[current_color][target_color].size = 
      neighb_list.size();
   if (offd_flag){
      MC_info->NeighbSend[current_color][target_color].ranks =
         (int *)calloc(MC_info->NeighbSend[current_color][target_color].size, 
                                                                sizeof(int));
      for (int q = 0; q < 
              MC_info->NeighbSend[current_color][target_color].size; q++){
         neighb_rank = neighb_list.front();

         neighb_list.pop_front();
         MC_info->NeighbSend[current_color][target_color].ranks[q] = 
                                                        neighb_rank;

         bound_ind_list[neighb_rank].sort();
         (*unique_bound_ind_list)[neighb_rank].assign(
            bound_ind_list[neighb_rank].begin(),
            bound_ind_list[neighb_rank].end());
         (*unique_bound_ind_list)[neighb_rank].unique();
      }
      MC_info->B[current_color][target_color] =
         (CSC *)malloc(MC_info->NeighbSend[current_color][target_color].size * 
                                                                 sizeof(CSC));
      for (int q = 0; q < 
           MC_info->NeighbSend[current_color][target_color].size; q++){
         neighb_rank = 
            MC_info->NeighbSend[current_color][target_color].ranks[q];
         MC_info->B[current_color][target_color][q].n = 
            MC_info->D[current_color].n;
         List_to_Block(&(MC_info->B[current_color][target_color][q]), 
                       bound_ind_list_2d[neighb_rank],
                       bound_elem_list_2d[neighb_rank]);
      }
   }
}


void DMEM_RMA_CM_SetupBoundWin_SymMat(DMEM_MulticolorInfo *MC_info,
                                      DMEM_ProcInfo Proc,
                                      vector<vector<vector<list<int>>>> 
                                         unique_bound_ind_list)
{
   int row, ind, shift_row;
   int neighb_rank, neighb_disp;
   int k, h, s;
   int map;
   int *sendbuff, *recvbuff;
   int *sendcnts, *recvcnts;
   int *sdispls, *rdispls;
   double elem;

   
   vector<list<int>>
      recv_disp_list(MC_info->P_glob.nparts);
   MC_info->Rma =
      (DMEM_RmaColorInfo *)malloc(MC_info->P_glob.nparts * sizeof(DMEM_RmaColorInfo));
   for (int c = 0; c < MC_info->P_glob.nparts; c++) MC_info->Rma[c].recv_count = 0;
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            neighb_rank = MC_info->NeighbSend[c][cc].ranks[q];
            MC_info->Bound[c][cc].num_send_points[q] = 
               unique_bound_ind_list[c][cc][neighb_rank].size();
            MC_info->Bound[c][cc].row_min[q] =
               unique_bound_ind_list[c][cc][neighb_rank].front();
            
            MC_info->Bound[c][cc].row_max[q] =
               unique_bound_ind_list[c][cc][neighb_rank].back();
            MC_info->Bound[c][cc].row_range_minmax[q] =
               (MC_info->Bound[c][cc].row_max[q] - 
                  MC_info->Bound[c][cc].row_min[q] + 1);
            MC_info->Bound[c][cc].row_range_disp[q] =
               (MC_info->Bound[c][cc].row_max[q] - 
                   MC_info->P_loc[cc].disp[neighb_rank] + 1);
            MC_info->Bound[c][cc].Rma.send_count[q] = 
               MC_info->Bound[c][cc].num_send_points[q];
            MC_info->Bound[c][cc].points[q] = 
               (int *)calloc(MC_info->Bound[c][cc].num_send_points[q], 
                                                         sizeof(int));
            MC_info->Bound[c][cc].Rma.send[q] =
               (double *)calloc(MC_info->Bound[c][cc].Rma.send_count[q], 
                                                        sizeof(double));
            for (int i = 0; i < 
                    MC_info->Bound[c][cc].num_send_points[q]; i++){
               MC_info->Bound[c][cc].points[q][i] = 
                  unique_bound_ind_list[c][cc][neighb_rank].front();
               unique_bound_ind_list[c][cc][neighb_rank].pop_front();
            }
            
            MC_info->Bound[c][cc].send_points[q] =
               (int *)calloc(MC_info->Bound[c][cc].row_range_minmax[q], 
                                                          sizeof(int));
            for (int i = 0; i < 
                    MC_info->Bound[c][cc].row_range_minmax[q]; i++)
               MC_info->Bound[c][cc].send_points[q][i] = -1;
            for (int i = 0; i < 
                 MC_info->Bound[c][cc].num_send_points[q]; i++){
               ind = 
                  MC_info->Bound[c][cc].points[q][i] - 
                     MC_info->Bound[c][cc].row_min[q];
               MC_info->Bound[c][cc].send_points[q][ind] = i;
               MC_info->Bound[c][cc].points[q][i] -= 
                  (MC_info->P_glob.disp[cc] + 
                     MC_info->P_loc[cc].disp[neighb_rank]);
               neighb_disp = 0;
               for (int ccc = 0; ccc < cc; ccc++)
                  neighb_disp += MC_info->P_loc[ccc].part[neighb_rank];
               MC_info->Bound[c][cc].points[q][i] += neighb_disp;
            }
         }

         sendbuff = (int *)calloc(Proc.world_size, sizeof(int));
         recvbuff = (int *)calloc(Proc.world_size, sizeof(int));
         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            neighb_rank = MC_info->NeighbSend[c][cc].ranks[q];
            sendbuff[neighb_rank] = 
               MC_info->Bound[c][cc].num_send_points[q];
         }
         MPI_Alltoall(sendbuff,
                      1,
                      MPI_INT,
                      recvbuff,
                      1,
                      MPI_INT,
                      MPI_COMM_WORLD); 
         MC_info->NeighbRecv[c][cc].size = 0;
         for (int p = 0; p < Proc.world_size; p++){
            if (recvbuff[p] > 0){
               MC_info->NeighbRecv[c][cc].size++;
            }
         }
         MC_info->NeighbRecv[c][cc].ranks = 
            (int *)calloc(MC_info->NeighbRecv[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].num_recv_points =
            (int *)calloc(MC_info->NeighbRecv[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].Rma.recv_disp =
            (int *)calloc(MC_info->NeighbRecv[c][cc].size, sizeof(int));
         MC_info->Bound[c][cc].recv_points =
            (int **)calloc(MC_info->NeighbRecv[c][cc].size,
                                            sizeof(int *));
         s = 0;
         for (int p = 0; p < Proc.world_size; p++){
            if (recvbuff[p] > 0){
               MC_info->NeighbRecv[c][cc].ranks[s] = p;
               s++;
            }
            sendbuff[p] = 0;
         }
         for (int q = 0; q < MC_info->NeighbRecv[c][cc].size; q++){
            neighb_rank = MC_info->NeighbRecv[c][cc].ranks[q];
            MC_info->Bound[c][cc].num_recv_points[q] = 
                                recvbuff[neighb_rank];
            recv_disp_list[c].push_back(MC_info->Rma[c].recv_count);
            sendbuff[neighb_rank] = MC_info->Rma[c].recv_count;
            MC_info->Rma[c].recv_count += 
               MC_info->Bound[c][cc].num_recv_points[q];
            MC_info->Bound[c][cc].recv_points[q] =
               (int *)calloc(MC_info->Bound[c][cc].num_recv_points[q],
                                                         sizeof(int));
         }
         
         for (int p = 0; p < Proc.world_size; p++){
            recvbuff[p] = 0;
         }
         MPI_Alltoall(sendbuff,
                      1,
                      MPI_INT,
                      recvbuff,
                      1,
                      MPI_INT,
                      MPI_COMM_WORLD);
         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            neighb_rank = MC_info->NeighbSend[c][cc].ranks[q];
            MC_info->Bound[c][cc].Rma.send_disp[q] =
               recvbuff[neighb_rank];
         }
         free(sendbuff);
         free(recvbuff);

         sendcnts = (int *)calloc(Proc.world_size, sizeof(int));
         recvcnts = (int *)calloc(Proc.world_size, sizeof(int));
         sdispls  = (int *)calloc((Proc.world_size+1), sizeof(int));
         rdispls  = (int *)calloc((Proc.world_size+1), sizeof(int));

         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            neighb_rank = MC_info->NeighbSend[c][cc].ranks[q];
            sendcnts[neighb_rank] =
               MC_info->Bound[c][cc].num_send_points[q];
         }
         for (int q = 0; q < MC_info->NeighbRecv[c][cc].size; q++){
            neighb_rank = MC_info->NeighbRecv[c][cc].ranks[q];
            recvcnts[neighb_rank] =
               MC_info->Bound[c][cc].num_recv_points[q];
         }
         for (int p = 0; p < Proc.world_size; p++){
            sdispls[p+1] = sdispls[p] + sendcnts[p];
            rdispls[p+1] = rdispls[p] + recvcnts[p];
         }

         sendbuff = (int *)calloc(sdispls[Proc.world_size], sizeof(int));
         recvbuff = (int *)calloc(rdispls[Proc.world_size], sizeof(int));

         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            neighb_rank = MC_info->NeighbSend[c][cc].ranks[q];
            for (int i = 0; i < 
                 MC_info->Bound[c][cc].num_send_points[q]; i++){
               sendbuff[sdispls[neighb_rank]+i] 
                  = MC_info->Bound[c][cc].points[q][i];
            }
         }
         MPI_Alltoallv(sendbuff,
                       sendcnts,
                       sdispls,
                       MPI_INT,
                       recvbuff,
                       recvcnts,
                       rdispls,
                       MPI_INT,
                       MPI_COMM_WORLD);
         for (int q = 0; q < MC_info->NeighbRecv[c][cc].size; q++){
            neighb_rank = MC_info->NeighbRecv[c][cc].ranks[q];
            for (int i = 0; i < 
                 MC_info->Bound[c][cc].num_recv_points[q]; i++){
               MC_info->Bound[c][cc].recv_points[q][i] = 
                  recvbuff[rdispls[neighb_rank]+i];
            }
         }

         free(sendbuff);
         free(recvbuff);
         free(sendcnts);
         free(recvcnts);
         free(sdispls);
         free(rdispls);
      }
   }
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
       MC_info->Rma[c].recv_disp_size = recv_disp_list[c].size();
       MC_info->Rma[c].recv_disp
          = (int *)calloc(MC_info->Rma[c].recv_disp_size, sizeof(int));
       for (int i = 0; i < MC_info->Rma[c].recv_disp_size; i++){
          MC_info->Rma[c].recv_disp[i] = recv_disp_list[c].front(); 
          recv_disp_list[c].pop_front();
       }
       MPI_Win_allocate(sizeof(double) * MC_info->Rma[c].recv_count,
                        sizeof(double),
                        MPI_INFO_NULL,
                        MPI_COMM_WORLD,
                        &(MC_info->Rma[c].recv),
                        &(MC_info->Rma[c].win));
       for (int i = 0; i < MC_info->Rma[c].recv_count; i++){
          MC_info->Rma[c].recv[i] = 0;
       }
   }
   MPI_Barrier(MPI_COMM_WORLD);

  // if (Proc.rank == 0){
  //    printf("rank %d:\n", Proc.rank);
  //    for (int c = 0; c < MC_info->P_glob.nparts; c++){
  //       for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
  //          printf("\t(%d, %d, %d):\n\t\t", c, cc, MC_info->P_loc[c].part[Proc.rank]);
  //          for (int q = 0; q < MC_info->NeighbSend[cc][c].size; q++){
  //             neighb_rank = MC_info->NeighbSend[cc][c].ranks[q];
  //             for (int i = 0; i < MC_info->Bound[cc][c].num_send_points[q]; i++){
  //                printf("{%d, %d}, ", neighb_rank, MC_info->Bound[cc][c].send_points[q][i]); 
  //             }
  //          }
  //          printf("\n");
  //       }
  //    }
  //    printf("\n");
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (Proc.rank == 1){
  //    printf("rank %d:\n", Proc.rank);
  //    for (int c = 0; c < MC_info->P_glob.nparts; c++){
  //       for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
  //          printf("\t(%d, %d, %d):\n\t\t", c, cc, MC_info->P_loc[c].part[Proc.rank]);
  //          for (int q = 0; q < MC_info->NeighbRecv[cc][c].size; q++){
  //             neighb_rank = MC_info->NeighbRecv[cc][c].ranks[q];
  //             for (int i = 0; i < MC_info->Bound[cc][c].num_send_points[q]; i++){
  //                printf("{%d, %d}, ", neighb_rank, MC_info->Bound[cc][c].send_points[q][i]);
  //             }
  //          }
  //          printf("\n");
  //       }
  //    }
  //    printf("\n");
  // }
  // MPI_Barrier(MPI_COMM_WORLD);


   MPI_Group world_group;
   int *temp_ranks = (int *)calloc(Proc.world_size, sizeof(int));
   for (int p = 0; p < Proc.world_size; p++)
      temp_ranks[p] = p;
   MPI_Comm_group(MPI_COMM_WORLD, &world_group);
   MPI_Group_incl(world_group, Proc.world_size,
                  temp_ranks, &(MC_info->mpi_group));
   MPI_Barrier(MPI_COMM_WORLD);
   free(temp_ranks);

}

void DMEM_RowPart(OrderInfo *P,
                  int world_size,
                  int extra_rank,
                  int n)
{
   int diff, extra, world, diff_prev, world_prev, extra_prev;

   if (floor(n / (double)world_size) != (n / (double)world_size)){
      world = (n - (n % world_size)) / world_size;
      extra = n - world * (world_size - 1);
      diff_prev = fabs(world - extra);
      while (1){
         world_prev = world;
         extra_prev = extra;
         if (world > extra){
            world--;
            extra += world_size-1;
         }
         else {
            world++;
            extra -= world_size-1;
         }
         diff = fabs(world - extra);
         if ((diff >= diff_prev) || (extra <= 0) || (world <= 0)){
            if ((extra <= 0) || (world <= 0)){
               world = world_prev;
               extra = extra_prev;
            }
            break;
         }
         diff_prev = diff;
      }
   }
   else{
      extra = world = n / world_size;
   }

   P->disp = (int *)malloc((world_size+1) * sizeof(int));
   P->part = (int *)malloc(world_size * sizeof(int));

   for (int i = 0; i < world_size; i++){
      P->disp[i] = i*world;
      if (i != extra_rank) P->part[i] = world;
      else P->part[i] = extra;
   }
   P->disp[world_size] = (world_size-1)*world + extra;
}

void DMEM_ZeroProc(DMEM_ProcInfo *Proc)
{
   for (int i = 0; i < Proc->Neighb.size; i++){
      if (solver_flag == DMEM_SOLVER_DS){
         Proc->Res.neighb_explicit_flag[i] = 0; 
         Proc->Res.neighb_norm_estim[i] = 0;
         Proc->Res.norm_estim_squ[i] = 0;
         Proc->Res.my_norm_estim_squ[i] = 0;
         Proc->Res.my_norm_estim_squ_prev[i] = 0;
         for (int j = 0; j < Proc->Bound.Rma.send_count[i]; j++){
            Proc->Bound.Rma.send_prev[i][j] = 0;
         }
         for (int j = 0; j < Proc->Bound.num_send_points[i]; j++){
            Proc->Bound.res[i][j] = 0;
         }
      }
      else if (solver_flag == DMEM_SOLVER_PS){
         Proc->Res.Rma.recv[i] = 0;
         Proc->Res.Rma.recv_prev[i] = 0;
      }
      for (int j = 0; j < Proc->Bound.Rma.send_count[i]; j++){
         Proc->Bound.Rma.send[i][j] = 0;
      }
   }
   for (int i = 0; i < Proc->Neighb.size+1; i++){
      if (solver_flag == DMEM_SOLVER_DS){
         Proc->Res.norm[i] = 0;
         Proc->Res.norm_prev[i] = 0;
      }
      else if (solver_flag == DMEM_SOLVER_PS){
         Proc->Res.norm[i] = 0;
         Proc->Res.norm_prev[i] = 0;
      }
   }
   for (int i = 0; i < Proc->Bound.Rma.recv_count; i++){
      Proc->Bound.Rma.recv[i] = 0;
      Proc->Bound.Rma.recv_prev[i] = 0;
   }
   if (implement_flag == DMEM_POS){
      for (int i = 0; i < Proc->Conv.win_size; i++){
         Proc->Conv.recv[i] = 0;
      }
      Proc->Conv.my_flag = 0;
      Proc->Conv.phase = 1;
   }
   if (solver_flag == DMEM_SOLVER_PS){
       Proc->Res.norm_send_prev[0] = 0;
   }
}

