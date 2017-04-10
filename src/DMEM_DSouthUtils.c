#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "DMEM_SolveUtils.h"
#include "DMEM_DSouthUtils.h"
#include "SEQ_Solve.h"

int DMEM_DSOUTH_RecvDiff(DMEM_ProcInfo *Proc, int q, int k)
{
   double diff;
   double x_elem, y_elem;
   int n = Proc->Bound.num_send_points[q] + 2;
   int return_val = 0;

   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Bound.Rma.win);
   for (int i = 0; i < n; i++){
      x_elem = Proc->Bound.Rma.recv[k+i];
      y_elem = Proc->Bound.Rma.recv_prev[k+i];

      if (!AlmostEqual2sComplement_DBL(x_elem, y_elem, 1)){
         return_val = 1;
         break;
      } 
   }
   MPI_Win_unlock(Proc->rank, Proc->Bound.Rma.win);
   return return_val;
}


void DMEM_DSOUTH_OffDiagSweep(DMEM_MatrixInfo Mat,
                              DMEM_ProcInfo *Proc,
                              SolveData *Out_data,
                              SolveVars *Vars)
{
   int ind, row, map, shift_row, k, neighb_rank, s;
   double elem, relax_val, start;

   start = MPI_Wtime();
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int i = 0; i < Mat.n; i++){
         k = Mat.B[q].j_ptr[i];
         for (int j = 0; j < Mat.B[q].j_ptr[i+1]-Mat.B[q].j_ptr[i]; j++) {
            ind = k+j;
            elem = Mat.B[q].a[ind];
            row = Mat.B[q].i[ind];
            shift_row = row - Proc->Bound.row_min[q];
            map = Proc->Bound.send_points[q][shift_row];
            relax_val = Vars->du[i] * elem;

            Proc->Bound.Rma.send[q][map] += relax_val;
            Proc->Bound.res[q][map] -= relax_val;
         }
      }
   }
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->wtime_update_bound[0] += (MPI_Wtime() - start);
}

void DMEM_DSOUTH_UpdateResEstim_Norm2(DMEM_MatrixInfo Mat,
                                      DMEM_ProcInfo *Proc,
                                      SolveData *Out_data)
{
   double start;
   start = MPI_Wtime();
   for (int q = 0; q < Proc->Neighb.size; q++){
      Proc->Res.my_norm_estim_squ[q] =
            SumSquaredDouble(Proc->Bound.res[q], 
                             Proc->Bound.num_send_points[q]);
      Proc->Res.norm_estim_squ[q] +=
         (Proc->Res.my_norm_estim_squ[q] -
             Proc->Res.my_norm_estim_squ_prev[q]);
      Proc->Res.my_norm_estim_squ_prev[q] = Proc->Res.my_norm_estim_squ[q];
      Proc->Res.norm[q+1] = sqrt(Proc->Res.norm_estim_squ[q]);
   }
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->wtime_res_estim[0] += (MPI_Wtime() - start);
}

void DMEM_RMA_DSOUTH_GatherBoundRes(DMEM_MatrixInfo Mat,
                                    DMEM_ProcInfo *Proc,
                                    DMEM_SolveVars Vars,
                                    SolveData *Out_data,
                                    int q)
{
   int k = 0, kk, point;
   double r_norm_squ, r_bound_norm_squ, start;

   start = MPI_Wtime();
   for (int i = 0; i < Proc->Bound.num_recv_points[q]; i++){
      point = Proc->Bound.recv_points[q][i];
      kk = Proc->Bound.num_send_points[q]+k;
      /* Accum boundary residual components.  Because we are
       * using accumulate, subtract out previous residuals. */
      Proc->Bound.Rma.send[q][kk] = (Vars.Loc.r[point] -
                      Proc->Bound.Rma.send_prev[q][kk]);
      Proc->Bound.Rma.send_prev[q][kk] = Vars.Loc.r[point];
      k++;
   }
   k = Proc->Bound.Rma.send_count[q]-2;
  // r_bound_norm_squ = Proc->Bound.res_norm;//SumSquaredDouble(Vars.Loc.r, Mat.n);
   Proc->Bound.Rma.send[q][k] =
      (pow(Proc->Bound.res_norm,2) - Proc->Bound.Rma.send_prev[q][k]);
   Proc->Bound.Rma.send_prev[q][k] = pow(Proc->Bound.res_norm,2);
   k++;

   Proc->Bound.Rma.send[q][k] =
      (Proc->Res.norm[q+1] - Proc->Bound.Rma.send_prev[q][k]);
   Proc->Bound.Rma.send_prev[q][k] = Proc->Res.norm[q+1];
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
}

void DMEM_RMA_DSOUTH_RecvResBound(DMEM_MatrixInfo Mat,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveParams Params,
                                  int q,
                                  int k)
{
   double elem;

   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Bound.Rma.win);
   for (int j = 0; j < Proc->Bound.num_send_points[q]; j++){
      Proc->Bound.res[q][j] = Proc->Bound.Rma.recv[k];
      Proc->Bound.Rma.recv_prev[k] = Proc->Bound.Rma.recv[k];
      k++;
   }
   Proc->Res.norm_estim_squ[q] = Proc->Bound.Rma.recv[k];
   Proc->Bound.Rma.recv_prev[k] = Proc->Bound.Rma.recv[k];
   k++;

   Proc->Res.neighb_norm_estim[q] = Proc->Bound.Rma.recv[k];
   Proc->Bound.Rma.recv_prev[k] = Proc->Bound.Rma.recv[k];
   k++;
   MPI_Win_unlock(Proc->rank, Proc->Bound.Rma.win);

   Proc->Res.norm[q+1] = sqrt(Proc->Res.norm_estim_squ[q]);
   Proc->Res.my_norm_estim_squ_prev[q] =
      SumSquaredDouble(Proc->Bound.res[q], Proc->Bound.num_send_points[q]);
}

void DMEM_RMA_DSOUTH_RecvAllResBound(DMEM_MatrixInfo Mat,
                                     DMEM_ProcInfo *Proc,
                                     DMEM_SolveVars *Vars,
                                     SolveParams Params)
{
   int map;
   double elem, diff, recv_norm, sum;
   int k = 0;

   for (int q = 0; q < Proc->Neighb.size; q++){
      k += Proc->Bound.num_recv_points[q];
      DMEM_RMA_DSOUTH_RecvResBound(Mat,
                                   Proc,
                                   Vars,
                                   Params,
                                   q,
                                   k);

      k += (Proc->Bound.num_send_points[q] + 2);
   }
}

void DMEM_RMA_DSOUTH_RecvAllBound(DMEM_MatrixInfo Mat,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveParams Params)
{
   int map, recv_flag;
   double elem, recv_norm, recv_diff, recv_diff_norm;
   int k = 0;

   for (int q = 0; q < Proc->Neighb.size; q++){
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Bound.Rma.win);
      for (int j = 0; j < Proc->Bound.num_recv_points[q]; j++){
         map = Proc->Bound.recv_points[q][j];
         elem = Proc->Bound.Rma.recv[k];
         Vars->Loc.r[map] -= elem;
         Proc->Bound.Rma.recv[k] = 0;
         k++;
      }
      MPI_Win_unlock(Proc->rank, Proc->Bound.Rma.win);

      recv_flag = DMEM_DSOUTH_RecvDiff(Proc, q, k);
      if (recv_flag){
         DMEM_RMA_DSOUTH_RecvResBound(Mat,
                                      Proc,
                                      Vars,
                                      Params,
                                      q,
                                      k);
      }
      k += (Proc->Bound.num_send_points[q] + 2);
   }
}

void DMEM_RMA_DSOUTH_AccumExplicitRes(DMEM_ProcInfo *Proc,
                                      SolveData *Out_data,
                                      int q)
{
   int target = Proc->Neighb.ranks[q];
   int k = Proc->Bound.num_recv_points[q];
   int count = Proc->Bound.num_send_points[q] + 2;
   int target_disp = Proc->Bound.Rma.send_disp[q] + k;

   MPI_Accumulate(&Proc->Bound.Rma.send[q][k],
                  count,
                  MPI_DOUBLE,
                  target,
                  target_disp,
                  count,
                  MPI_DOUBLE,
                  MPI_SUM,
                  Proc->Bound.Rma.win);
}

void DMEM_POS_DSOUTH_AccumExplicitRes(DMEM_ProcInfo *Proc,
                                      SolveData *Out_data,
                                      int q)
{
   double start;
   int target = Proc->Neighb.ranks[q];

   start = MPI_Wtime();
   MPI_Win_lock(MPI_LOCK_SHARED, target, 0, Proc->Bound.Rma.win);
   DMEM_RMA_DSOUTH_AccumExplicitRes(Proc, Out_data, q);
   MPI_Win_unlock(target, Proc->Bound.Rma.win);
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
}

void DMEM_DSOUTH_ExplicitResUpdate(DMEM_MatrixInfo Mat,
                                   DMEM_ProcInfo *Proc,
                                   SolveData *Out_data,
                                   DMEM_SolveVars *Vars,
                                   SolveParams Params,
                                   VSLStreamStatePtr stream)
{
   double rand_num, prob, start;
   int count = 0;
   bool equal_flag;

   if (Proc->Res.explicit_flag){
      for (int q = 0; q < Proc->Neighb.size; q++){
         if (Proc->Res.neighb_explicit_flag[q]){
            equal_flag = 
               AlmostEqual2sComplement_DBL(Proc->Bound.res_norm,
                                           Proc->Res.neighb_norm_estim[q],
                                           1);
            if ((Proc->Bound.res_norm < Proc->Res.neighb_norm_estim[q]) &&
                (equal_flag == false)) {
               DMEM_RMA_DSOUTH_GatherBoundRes(Mat, 
                                              Proc, 
                                              *Vars, 
                                              Out_data,
                                              q);
               if (implement_flag == DMEM_POS){
                  DMEM_POS_DSOUTH_AccumExplicitRes(Proc, Out_data, q);
                 // DMEM_POS_AccumBound(Proc, Out_data, q);
               }
               else if(implement_flag == DMEM_SOS){
                  DMEM_RMA_DSOUTH_AccumExplicitRes(Proc, Out_data, q);
                 // DMEM_RMA_AccumBound(Proc, q);
               }
               Out_data->comm[0]++;
               Out_data->comm_res[0]++;
               Out_data->comm_scaled[0] += (1.0/(double)Proc->world_size);
               count++;
            }
            Proc->Res.neighb_explicit_flag[q] = 0;
         }
      }
      Proc->Res.explicit_flag = 0;
   }
}

void DMEM_DSOUTH_SOS_ExplicitResUpdate(DMEM_MatrixInfo Mat,
                                       DMEM_ProcInfo *Proc,
                                       SolveData *Out_data,
                                       DMEM_SolveVars *Vars,
                                       SolveParams Params,
                                       VSLStreamStatePtr stream)
{
   double start;
   start = MPI_Wtime();
   MPI_Win_post(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
   MPI_Win_start(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);

   DMEM_DSOUTH_ExplicitResUpdate(Mat, Proc, Out_data, Vars, Params, stream);

   MPI_Win_complete(Proc->Bound.Rma.win);
   MPI_Win_wait(Proc->Bound.Rma.win);
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
}
