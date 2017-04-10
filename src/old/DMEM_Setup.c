#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Laplace.h"
#include "Multicolor.h"
#include "FileUtils.h"
#include "MatrixUtils.h"
#include "Misc.h"

void DMEM_DistributeMatrix(CSC A,
                           DMEM_MatrixInfo *Mat,
                           DMEM_ProcInfo *Proc);
void DMEM_SetupBoundary(DMEM_MatrixInfo *Mat, 
                        DMEM_ProcInfo *Proc);
void DMEM_RMA_SetupUniqueBoundWin_NonSymMat(DMEM_ProcInfo *Proc,
                                            DMEM_MatrixInfo *Mat,
                                            std::vector<std::list<int>>
                                               unique_bound_ind_list);
void DMEM_RMA_SetupBoundWin_NonSymMat(DMEM_ProcInfo *Proc,
                                      DMEM_MatrixInfo *Mat,
                                      std::vector<std::list<int>>
                                         bound_ind_list);
void DMEM_RMA_SetupBoundWin_SymMat(DMEM_ProcInfo *Proc,
                                   DMEM_MatrixInfo *Mat,
                                   std::vector<std::list<int>>
                                      bound_ind_list);

void DMEM_Setup(FILE *in_file,
                DMEM_MatrixInfo *Mat,
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

   idx_t nparts;
   if (color_flag){
      nparts = 2;
   }
   else {
      nparts = (idx_t)(Proc->world_size + imbalance);
      Mat->P.nparts = Proc->world_size;
   }

   idx_t options[METIS_NOPTIONS];
   char buffer[256];
   FILE *temp_file;
   double start, stop;

   if ((Proc->rank == 0) && (!format_out_flag))
      printf("\n|********** INITIALIZING PROBLEM **********|\n");
   start = MPI_Wtime();
   if (mat_file_flag){
      if (Proc->rank == 0){
         if (!format_out_flag){
            printf("\nReading matrix from file.\n");
         }
         ReadBinary_fread_metis(in_file, &G, &T);
      }
   }
   else {
      if (Proc->rank == 0){
         if (!format_out_flag){
            printf("\nUsing 5-pt centered difference discretization of the Laplace PDE.\n");
         }
         Laplace2D_FD5pt_metis(&G, &T, m, w); 
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   stop = MPI_Wtime() - start;
   if ((Proc->rank == 0) && !format_out_flag){
      printf("\nMatrix loaded, time = %es, nnz = %d, n = %d.\n", 
             stop, G.nnz, G.n);
   }

   METIS_SetDefaultOptions(options);
   options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

   A.n = G.n;
   A.nnz = G.nnz;
   if (nparts > 1){
      if (Proc->rank == 0){
         if (color_flag){
           // if (Proc->rank == 0){
           //    Metis_to_CSC(&A, G);
           //    if (!format_out_flag) printf("\nREORDERING USING MULTICOLORING\n");
           //    Multicolor(Mat->A, &(Mat->P));
           // }
         }
         else {
            idx_t *perm = (idx_t *)calloc(A.n, sizeof(idx_t));
            start = MPI_Wtime();
            flag =  METIS_PartGraphKway(&(G.n), &ncon, G.xadj, G.adjncy, NULL, 
                                        NULL, NULL, &nparts, NULL, NULL, 
                                        options, &objval, perm);
            stop = MPI_Wtime() - start;
            if (flag != METIS_OK)
               printf("****WARNING****: METIS returned error with code %d.\n", flag);
            if (!format_out_flag) 
               printf("\nPartitions determined using METIS, time = %es.\n", stop);
           // Free_metis(&G); 
            Mat->P.perm = (int *)malloc(A.n * sizeof(int));
            Mat->P.part = (int *)calloc(Mat->P.nparts, sizeof(int));
            for (int i = 0; i < A.n; i++){
               Mat->P.perm[i] = perm[i];
               Mat->P.part[Mat->P.perm[i]]++;
            }
            Mat->P.disp = (int *)malloc((Mat->P.nparts+1) * sizeof(int));
            Mat->P.disp[0] = 0;
            for (int i = 0; i < Mat->P.nparts; i++){
               Mat->P.disp[i+1] = Mat->P.disp[i] + Mat->P.part[i];
            }
            start = MPI_Wtime();
            ReorderTriplet_csc(T, &A, &(Mat->P));
            stop = MPI_Wtime() - start;
            free(T.a);
            free(T.i);
            free(T.j);
            if (!format_out_flag)
               printf("\nMatrix reordered, time = %es.\n", stop);
           // Write_csc(A, 1);
         }
      }
     // if (color_flag){
     //    if (!format_out_flag){
     //        printf("colors = %d\n"
     //               "partition: mean = %.2f, max %d, min %d\n",
     //               Mat->P.nparts, MeanInt(Mat->P.part, Mat->P.nparts),
     //               MaxInt(Mat->P.part,nparts), MinInt(Mat->P.part,nparts));
     //        printf("%d & %d & %d & %d & %d\\\\",
     //               nnz, n, Mat->P.nparts, MaxInt(Mat->P.part,nparts),
     //               MinInt(Mat->P.part,nparts));
     //    }
     // }
      start = MPI_Wtime();
      DMEM_DistributeMatrix(A, Mat, Proc);
      MPI_Barrier(MPI_COMM_WORLD);
      stop = MPI_Wtime() - start;
      if ((Proc->rank == 0) && !format_out_flag) 
         printf("\nMatrix distributed, time = %es.\n", stop);
      start = MPI_Wtime();
      DMEM_SetupBoundary(Mat, Proc);
      MPI_Barrier(MPI_COMM_WORLD);
      stop = MPI_Wtime() - start;
      if ((Proc->rank == 0) && !format_out_flag)
         printf("\nBoundary point information initialized, time = %es.\n", stop);
     // DMEM_WriteBlocks_csc(*Mat, *Proc, 1);
   }
   else{
      Mat->P.dispv = (int *)calloc(1, sizeof(int));
      Mat->disp = Mat->P.dispv[0] = 0;
      Mat->P.part = (int *)calloc(1, sizeof(int));
      Mat->n = Mat->n_glob = Mat->A.n = Mat->P.part[0] = A.n;
      Mat->nnz = Mat->A.nnz = A.nnz;
      Metis_to_CSC(&Mat->A, G);
      Free_metis(&G);
   }
  // Write_csc(Mat->A, 1);
   if ((Proc->rank == 0) && !format_out_flag) 
      printf("\n|******************************************|\n\n");
}

void DMEM_DistributeMatrix(CSC A,
                           DMEM_MatrixInfo *Mat, 
                           DMEM_ProcInfo *Proc)
{
   int disp;
   int j_ptr_disp;
   Mat->P.dispv = (int *)malloc(Mat->P.nparts * sizeof(int));
   int *j_ptr_extra = (int *)malloc(Mat->P.nparts * sizeof(int));
   if (Proc->rank == 0){
      Mat->n_glob = A.n;
      for (int i = 0; i < Mat->P.nparts; i++){ 
         Mat->P.dispv[i] = Mat->P.disp[i];
         j_ptr_extra[i] = A.j_ptr[Mat->P.disp[i+1]];
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
   MPI_Scatterv(A.diag, Mat->P.part, Mat->P.dispv, MPI_DOUBLE,
                Mat->A.diag, Mat->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   Mat->A.j_ptr = (int *)calloc(Mat->n+1, sizeof(int));
   MPI_Scatterv(A.j_ptr, 
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
            count_nnz += (A.j_ptr[j+1] - A.j_ptr[j]);
         }
         if (i > 0) dispv_nnz[i] = dispv_nnz[i-1] + count_nnz_prev;
         part_nnz[i] = count_nnz_prev = count_nnz;
        // printf("%d, %d, %d\n", count_nnz, dispv_nnz[i], Mat->P.part[i]);
      }
   }
   
  // if (Proc->rank == 0){
  //    int row, col, k;
  //    double elem;
  //    char buffer[100];
  //    strcpy(buffer, "metis_matrix_matlab.txt");
  //    FILE *out_file;
  //    remove(buffer);
  //    out_file = fopen(buffer, "a");
  //    for (int i = 0; i < A.n; i++){
  //       k = A.j_ptr[i];
  //       for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++){
  //          row = A.i[k+j];
  //          col = i;
  //          elem = A.a[k+j];
  //          fprintf(out_file, "%d %d %e\n",
  //                  row+1, col+1, elem);
  //       }
  //    }
  // }
   
   MPI_Bcast(part_nnz, Mat->P.nparts, MPI_INT, 0, MPI_COMM_WORLD);
   Mat->A.nnz = Mat->nnz = part_nnz[Proc->rank];
   Mat->A.i = (int *)calloc(Mat->A.nnz, sizeof(int));
   Mat->A.a = (double *)calloc(Mat->A.nnz, sizeof(double));
   MPI_Scatterv(A.i, part_nnz, dispv_nnz, MPI_INT,
                Mat->A.i, Mat->A.nnz, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Scatterv(A.a, part_nnz, dispv_nnz, MPI_DOUBLE,
                Mat->A.a, Mat->A.nnz, MPI_DOUBLE, 0,MPI_COMM_WORLD);

   free(part_nnz);
   free(dispv_nnz);
   free(j_ptr_extra);
}

void DMEM_SetupBoundary(DMEM_MatrixInfo *Mat,
                        DMEM_ProcInfo *Proc)
{
   int cum_part, next_part;
   int row_rank;
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   double elem;

   /* Sort columns of un-split matrix (required by PARDISO). */
   for (int i = 0; i < Mat->n; i++){
      QuicksortPair_int_dbl(Mat->A.i, 
                            Mat->A.a, 
                            Mat->A.j_ptr[i], 
                            Mat->A.j_ptr[i+1]-1);
   }

   std::vector<std::list<int>> diag_ind_list(Mat->n);
   std::vector<std::list<double>> diag_elem_list(Mat->n);
   std::vector<std::list<int>> unique_bound_ind_list(Proc->world_size);
   std::vector<std::list<int>> bound_ind_list(Proc->world_size);
   std::vector<std::vector<std::list<int>>> 
      bound_ind_list_2d(Proc->world_size, 
                        std::vector<std::list<int>>(Mat->n));
   std::vector<std::vector<std::list<double>>>
      bound_elem_list_2d(Proc->world_size,
                         std::vector<std::list<double>>(Mat->n));
   std::list<int> neighb_list;
   for (int i = 0; i < Mat->n; i++){
      for (int j = 0; j < Mat->A.j_ptr[i+1]-Mat->A.j_ptr[i]; j++){
         k = Mat->A.j_ptr[i]+j;
         row = Mat->A.i[k];
         elem = Mat->A.a[k];
         row_rank = 0;
         while (Mat->P.disp[row_rank] < Mat->n_glob){
            next_part = Mat->P.part[row_rank];
            if (row < (Mat->P.disp[row_rank] + next_part)) {
               if (row_rank == Proc->rank) {
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
      unique_bound_ind_list[neighb_rank].assign(
         bound_ind_list[neighb_rank].begin(), 
         bound_ind_list[neighb_rank].end());
      unique_bound_ind_list[neighb_rank].unique();
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

   Proc->Res.norm1 =
      (double *)calloc(Proc->Neighb.size+1, sizeof(double *));
   Proc->Bound.size =
       (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.row_min =
      (int *)calloc(Proc->Neighb.size, sizeof(int *));
   Proc->Bound.row_max =
      (int *)calloc(Proc->Neighb.size, sizeof(int *));
   Proc->Bound.Rma.send_map =
      (int **)calloc(Proc->Neighb.size, sizeof(int *));
   Proc->Bound.unique_size =
      (int *)calloc(Proc->world_size, sizeof(int));
   Proc->Bound.unique_range =
      (int *)calloc(Proc->world_size, sizeof(int));
   Proc->Bound.Rma.send_count =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.send_disp =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.recv_part =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.recv_disp =
      (int *)calloc(Proc->Neighb.size, sizeof(int));
   Proc->Bound.Rma.send =
      (double **)calloc(Proc->Neighb.size, sizeof(double *));
   Proc->Bound.Rma.recv_map =
      (int **)calloc(Proc->Neighb.size, sizeof(int *));

   if (solver_flag == DMEM_SOLVER_DS){
      Proc->Bound.Rma.send_prev =
      (double **)calloc(Proc->Neighb.size, sizeof(double *));
   }
    

   DMEM_RMA_SetupUniqueBoundWin_NonSymMat(Proc, 
                                          Mat, 
                                          unique_bound_ind_list);
  // DMEM_RMA_SetupBoundWin_SymMat(Proc, 
  //                               Mat, 
  //                               bound_ind_list);
}

void DMEM_RMA_SetupUniqueBoundWin_NonSymMat(DMEM_ProcInfo *Proc,
                                            DMEM_MatrixInfo *Mat,
                                            std::vector<std::list<int>>
                                               unique_bound_ind_list)
{
   int row, ind;
   int neighb_rank;
   int k;
   int map;
   double elem;
   MPI_Request req_send, req_recv;
   MPI_Status stat_send, stat_recv;

   int **unique_bound_ind =
         (int **)calloc(Proc->Neighb.size, sizeof(int *)); 
   Proc->Bound.Rma.recv_count = 0;
   for (int i = 0; i < Proc->Neighb.size; i++){
      neighb_rank = Proc->Neighb.ranks[i];
      Proc->Bound.unique_size[i] = unique_bound_ind_list[neighb_rank].size();
      Proc->Bound.Rma.send_count[i] = Proc->Bound.unique_size[i];
      if (solver_flag == DMEM_SOLVER_DS){ 
         Proc->Bound.Rma.send_prev[i] =
            (double *)calloc(Proc->Bound.Rma.send_count[i], sizeof(double));
         Proc->Bound.Rma.send_count[i] += 2;
      }
      Proc->Bound.Rma.send[i] =
         (double *)calloc(Proc->Bound.Rma.send_count[i], sizeof(double));
      Proc->Bound.row_min[i] = unique_bound_ind_list[neighb_rank].front();
      Proc->Bound.row_max[i] = unique_bound_ind_list[neighb_rank].back();
      Proc->Bound.unique_range[i] =
        // (Proc->Bound.row_max[i] - Proc->Bound.row_min[i] + 1);
         (Proc->Bound.row_max[i] - Mat->P.disp[neighb_rank] + 1);
      unique_bound_ind[i] = 
         (int *)calloc(Proc->Bound.unique_size[i], sizeof(int));
      for (int j = 0; j < Proc->Bound.unique_size[i]; j++){
         unique_bound_ind[i][j] = unique_bound_ind_list[neighb_rank].front();
         unique_bound_ind[i][j] -= Mat->P.disp[neighb_rank];//Proc->Bound.row_min[i];
         unique_bound_ind_list[neighb_rank].pop_front();
      }
      
      Proc->Bound.Rma.send_map[i] =
         (int *)calloc(Proc->Bound.unique_range[i], sizeof(int));
      for (int j = 0; j < Proc->Bound.unique_range[i]; j++)
         Proc->Bound.Rma.send_map[i][j] = -1;
      for (int j = 0; j < Proc->Bound.unique_size[i]; j++){
         ind = unique_bound_ind[i][j];
         Proc->Bound.Rma.send_map[i][ind] = j;
      }
      MPI_Isend(&(Proc->Bound.unique_size[i]), 
                1, 
                MPI_INT, 
                neighb_rank, 
                0, 
                MPI_COMM_WORLD,
                &req_send);
      MPI_Recv(&(Proc->Bound.Rma.recv_part[i]), 
               1, 
               MPI_INT, 
               neighb_rank, 
               0, 
               MPI_COMM_WORLD, 
               &stat_recv);
      MPI_Wait(&req_send, &stat_recv);

      Proc->Bound.Rma.recv_disp[i] = Proc->Bound.Rma.recv_count;
      Proc->Bound.Rma.recv_count += Proc->Bound.Rma.recv_part[i];
      if (solver_flag == DMEM_SOLVER_DS){
         Proc->Bound.Rma.recv_count += 2;
      }      
      Proc->Bound.Rma.recv_map[i] =
         (int *)calloc(Proc->Bound.Rma.recv_part[i], sizeof(int));

      neighb_rank = Proc->Neighb.ranks[i];
      MPI_Isend(unique_bound_ind[i], 
                Proc->Bound.unique_size[i],
                MPI_INT, 
                neighb_rank,
                0, 
                MPI_COMM_WORLD,
                &req_send);
      MPI_Recv(Proc->Bound.Rma.recv_map[i], 
               Proc->Bound.Rma.recv_part[i], 
               MPI_INT, 
               neighb_rank,
               0, 
               MPI_COMM_WORLD, 
               &stat_recv);
      MPI_Wait(&req_send, &stat_send);

      MPI_Isend(&(Proc->Bound.Rma.recv_disp[i]),
                1,
                MPI_INT,
                neighb_rank,
                0,
                MPI_COMM_WORLD,
                &req_send);
      MPI_Recv(&(Proc->Bound.Rma.send_disp[i]),
               1,
               MPI_INT,
               neighb_rank,
               0,
               MPI_COMM_WORLD,
               &stat_recv);
      MPI_Wait(&req_send, &stat_send);
   }

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
   
   MPI_Win_allocate(sizeof(int),
                    sizeof(int),
                    MPI_INFO_NULL,
                    MPI_COMM_WORLD,
                    &(Proc->Conv.recv),
                    &(Proc->Conv.win)); 
   Proc->Conv.flag = 0; 
}

//void DMEM_RowPart(ProcInfo *proc)
//{
//   int rank_p, num_p, diff, extra, world, diff_prev;
//   int row_flag, sum;
//   int n = A->n;
//
//   MPI_Comm_size(MPI_COMM_WORLD, &num_p);
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank_p);
//
//   A->disp = malloc((num_p+1) * sizeof(int));
//   A->disp_v = malloc(num_p * sizeof(int));
//   A->part = calloc(num_p, sizeof(int));
//
//   if (floor(n / (double)num_p) != (n / (double)num_p)){
//      world = (n - (n % num_p)) / num_p;
//      extra = n - world * (num_p - 1);
//      diff_prev = fabs(world - extra);
//      while (1){
//         if (world > extra){
//            world--;
//            extra += num_p-1;
//         }
//         else {
//            world++;
//            extra -= num_p-1;
//         }
//         diff = fabs(world - extra);
//         if ((diff >= diff_prev) || (extra <= 0) || (world <= 0)) break;
//         diff_prev = diff;
//      }
//   }
//   else{
//      extra = world = n / num_p;
//   }
//
//   for (int i = 0; i < num_p; i++){
//      A->disp[i] = A->disp_v[i] = i*world;
//      if (i != num_p-1) A->part[i] = world;
//      else A->part[i] = extra;
//   }
//   A->disp[num_p] = (num_p-1)*world + extra;
//   A->np = A->part[rank_p];
//   A->max_np = MaxInt(A->part, num_p);
//}

//void CSC_RowMap(CSC *A)
//{
//   int num_p = 1, k;
//   int n = A->n;
//   int *ascend = (int *)calloc(n, sizeof(int));
//   A->row_to_p =  (int *)calloc(n, sizeof(int));
//   for (int i = 0; i < n; i++) ascend[i] = i;
//
//   k = 0;
//   for (int i = 0; i < n; i++){
//      if (ascend[i] >= A->disp[k*A->loc_parts]) k++;
//      A->row_to_p[i] = k-1;
//   }
//
//   free(ascend);
//}

/* INCOMPLETE FUNCTION */
//void DMEM_RMA_SetupBoundWin_SymMat(DMEM_ProcInfo *Proc,
//                                   DMEM_MatrixInfo *Mat,
//                                   std::vector<std::list<int>> 
//                                      bound_ind_list)
//{
//   int row, ind, neighb_rank, k, map;
//   double elem;
//   MPI_Request req_send, req_recv;
//   MPI_Status stat_send, stat_recv;
//
//   Proc->Bound.Rma.recv_count = 0;
//   for (int i = 0; i < Proc->Neighb.size; i++){
//      neighb_rank = Proc->Neighb.ranks[i];
//      Proc->Bound.Rma.send_count[i] = Mat->B[i].nnz;
//      Proc->Bound.Rma.send[i] =
//         (double *)calloc(Proc->Bound.Rma.send_count[i], sizeof(double));
//      Proc->Bound.Rma.recv_part[i] = Mat->B[i].nnz;
//      Proc->Bound.Rma.recv_disp[i] = Proc->Bound.Rma.recv_count;
//      Proc->Bound.Rma.recv_count += Proc->Bound.Rma.recv_part[i];
//   }
//   for (int i = 0; i < Proc->Neighb.size; i++){
//      MPI_Isend(&(Proc->Bound.Rma.recv_disp[i]),
//                1,
//                MPI_INT,
//                neighb_rank,
//                0,
//                MPI_COMM_WORLD,
//                &req_send);
//      MPI_Recv(&(Proc->Bound.Rma.send_disp[i]),
//               1,
//               MPI_INT,
//               neighb_rank,
//               0,
//               MPI_COMM_WORLD,
//               &stat_recv);
//      MPI_Wait(&req_send, &stat_send); 
//   }
//   Proc->Bound.Rma.recv_map =
//      (int *)calloc(Proc->Bound.Rma.recv_count, sizeof(int));
//   k = 0;
//   for (int i = 0; i < Proc->Neighb.size; i++){
//      neighb_rank = Proc->Neighb.ranks[i];
//      for (int j = 0; j < Mat->B[i].nnz; j++){
//         Proc->Bound.Rma.recv_map[k] = (Mat->B[i].i[j] - Mat->P.disp[neighb_rank]);
//         k++;
//      }
//   }
//
//   MPI_Win_allocate(sizeof(double) * Proc->Bound.Rma.recv_count,
//                    sizeof(double),
//                    MPI_INFO_NULL,
//                    MPI_COMM_WORLD,
//                    &(Proc->Bound.Rma.recv),
//                    &(Proc->Bound.Rma.win));
//   Proc->Bound.Rma.recv_prev =
//      (double *)calloc(Proc->Bound.Rma.recv_count, sizeof(double));
//   for (int i = 0; i < Proc->Bound.Rma.recv_count; i++){
//      Proc->Bound.Rma.recv[i] = Proc->Bound.Rma.recv_prev[i] = 0;
//   }
//   
//   MPI_Win_allocate(sizeof(int),
//                    sizeof(int),
//                    MPI_INFO_NULL,
//                    MPI_COMM_WORLD,
//                    &(Proc->Conv.recv),
//                    &(Proc->Conv.win)); 
//   Proc->Conv.flag = 0; 
//}

void DMEM_SetupOutData(SolveData *Out_data)
{
   Out_data->wtime_tot = (double *)calloc(1, sizeof(double));
   Out_data->comm_scaled = (double *)calloc(1, sizeof(double));
   Out_data->relax_scaled = (double *)calloc(1, sizeof(double));
   Out_data->sweep =
      (unsigned long long *)calloc(1, sizeof(unsigned long long));
   Out_data->relax =
      (unsigned long long *)calloc(1, sizeof(unsigned long long));
   Out_data->comm =
      (unsigned long long *)calloc(1, sizeof(unsigned long long));
   Out_data->sweep[0] = 0;
   Out_data->relax[0] = 0;
   Out_data->comm[0] = 0;
}
