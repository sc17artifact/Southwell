#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "SEQ_Solve.h"

void DMEM_UpdateResEstim_Norm1(DMEM_MatrixInfo Mat, 
                               DMEM_ProcInfo *Proc)
{
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int i = 0; i < Proc->Bound.unique_size[q]; i++){
         Proc->Res.norm1[q+1] +=
            (Proc->Bound.Rma.send[q][i] -
               Proc->Bound.Rma.send_prev[q][i]);
      }
   }
}

void RMA_Recv(DMEM_MatrixInfo Mat,
              DMEM_ProcInfo *Proc,
              double *r)
{
   int map;
   double elem;

   int k = 0;
   for (int i = 0; i < Proc->Neighb.size; i++){
      for (int j = 0; j < Proc->Bound.Rma.recv_part[i]; j++){
         map = Proc->Bound.Rma.recv_map[i][j];
         #pragma omp atomic read
         elem = Proc->Bound.Rma.recv[k];
         r[map] -= elem;
         #pragma omp atomic write
         Proc->Bound.Rma.recv[k] = 0;
         k++;
      }
   }
}

void RMA_Recv_DistrSouthwell(DMEM_MatrixInfo Mat,
                             DMEM_ProcInfo *Proc,
                             double *r)
{
   int map;
   double elem;

   int k = 0;
   for (int i = 0; i < Proc->Neighb.size; i++){
      for (int j = 0; j < Proc->Bound.Rma.recv_part[i]; j++){
         map = Proc->Bound.Rma.recv_map[i][j];
         #pragma omp atomic read
         elem = Proc->Bound.Rma.recv[k];
         r[map] -= elem;
         #pragma omp atomic write
         Proc->Bound.Rma.recv[k] = 0;
         k++;
      }
      #pragma omp atomic read
      elem = Proc->Bound.Rma.recv[k];
      if (elem > 0){
         Proc->Res.norm1[i+1] = elem;
         #pragma omp atomic write
         Proc->Bound.Rma.recv[k] = 0;
      }
      k++;
   
      #pragma omp atomic read
      elem = Proc->Bound.Rma.recv[k];
      if (elem > 0){
         Proc->Res.norm1_estim[i] = elem;
         #pragma omp atomic write
         Proc->Bound.Rma.recv[k] = 0;
      }
      k++;
   }
}

void POS_SendAll(DMEM_ProcInfo *Proc)
{
   int target;
   for (int q = 0; q < Proc->Neighb.size; q++){
      target = Proc->Neighb.ranks[q];
      MPI_Win_lock(MPI_LOCK_SHARED, target, 0, Proc->Bound.Rma.win);
      MPI_Accumulate(Proc->Bound.Rma.send[q],
                     Proc->Bound.Rma.send_count[q],
                     MPI_DOUBLE,
                     target,
                     Proc->Bound.Rma.send_disp[q],
                     Proc->Bound.Rma.send_count[q],
                     MPI_DOUBLE,
                     MPI_SUM,
                     Proc->Bound.Rma.win);
      MPI_Win_unlock(target, Proc->Bound.Rma.win);
      for (int i = 0; i < Proc->Bound.Rma.send_count[q]; i++){
         Proc->Bound.Rma.send[q][i] = 0;
      }
   }
}

void POS_Send(DMEM_ProcInfo *Proc,
              int q)
{
   int target = Proc->Neighb.ranks[q];
   MPI_Win_lock(MPI_LOCK_SHARED, target, 0, Proc->Bound.Rma.win);
   MPI_Accumulate(Proc->Bound.Rma.send[q], 
                  Proc->Bound.Rma.send_count[q], 
                  MPI_DOUBLE,
                  target, 
                  Proc->Bound.Rma.send_disp[q], 
                  Proc->Bound.Rma.send_count[q], 
                  MPI_DOUBLE, 
                  MPI_SUM, 
                  Proc->Bound.Rma.win);
   MPI_Win_unlock(target, Proc->Bound.Rma.win);
   for (int i = 0; i < Proc->Bound.Rma.send_count[q]; i++){
      Proc->Bound.Rma.send[q][i] = 0;
   }
}

void DiagSweep_GS(DMEM_MatrixInfo Mat,
                  DMEM_ProcInfo Proc,
                  SolveData *Out_data,
                  SolveVars *Vars)
{
   int ind, row, shift_row;
   double u_prev, elem;
   double relax = 0;

   for (int i = 0; i < Mat.n; i++){
      u_prev = Vars->u[i];
      Vars->u[i] += (Vars->r[i] / Mat.D.diag[i]);
      Vars->du[i] = Vars->u[i] - u_prev;
      for (int j = 0; j < Mat.D.j_ptr[i+1]-Mat.D.j_ptr[i]; j++) {
         ind = Mat.D.j_ptr[i]+j;
         row = Mat.D.i[ind];
         elem = Mat.D.a[ind];
         shift_row = row - Mat.disp;
         Vars->r[shift_row] -= (Vars->du[i] * elem);
      }
      relax++;
      Out_data->relax[0]++;
   }
   Out_data->relax_scaled[0] += (relax / (double)Mat.n_glob);
}

void OffDiagSweep_GS(DMEM_MatrixInfo Mat,
                     DMEM_ProcInfo *Proc,
                     SolveData *Out_data,
                     SolveVars *Vars)
{
   int ind, row, map, shift_row, k, neighb_rank, s;
   double elem;
   double relax = 0;

   for (int q = 0; q < Proc->Neighb.size; q++){ 
      neighb_rank = Proc->Neighb.ranks[q];
      for (int i = 0; i < Mat.n; i++){
         k = Mat.B[q].j_ptr[i];
         for (int j = 0; j < Mat.B[q].j_ptr[i+1]-Mat.B[q].j_ptr[i]; j++) {
            ind = k+j;
            elem = Mat.B[q].a[ind];
            row = Mat.B[q].i[ind];
            shift_row = row - Mat.P.disp[neighb_rank];
            map = Proc->Bound.Rma.send_map[q][shift_row];
            Proc->Bound.Rma.send[q][map] += (Vars->du[i] * elem);
         }
         relax++;
         Out_data->relax[0]++;
      }
      Out_data->relax[0] += (relax / (double)Mat.n_glob);
   }
}

int POS_CheckConverge(double norm, 
                      SolveData Out_data, 
                      SolveParams Params,
                      DMEM_ProcInfo *Proc)
{
   int one = 1;
   if (((Out_data.sweep[0] >= Params.sweep_max) || (norm < Params.tol)) &&
       (!(Proc->Conv.flag))){
      for (int i = 0; i < Proc->world_size; i++){
         MPI_Win_lock(MPI_LOCK_SHARED, i, 0, Proc->Conv.win);
         MPI_Accumulate(&one, 1, MPI_INT, i, 0, 1, 
                        MPI_INT, MPI_SUM, Proc->Conv.win);
         MPI_Win_unlock(i, Proc->Conv.win);
      }
      Proc->Conv.flag = 1;
   }
   return Proc->Conv.recv[0];
}
