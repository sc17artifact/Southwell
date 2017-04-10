#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "SEQ_Solve.h"

extern int ds_fail, ds_pass;
extern int loc_solver_flag;

void DMEM_RMA_RecvBound(DMEM_MatrixInfo Mat,
                        DMEM_ProcInfo *Proc,
                        double *r)
{
   int map;
   double elem;

   int k = 0;
   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Bound.Rma.win);
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int j = 0; j < Proc->Bound.num_recv_points[q]; j++){
         map = Proc->Bound.recv_points[q][j];
         elem = Proc->Bound.Rma.recv[k];
         r[map] -= elem;
         Proc->Bound.Rma.recv[k] = 0;
         k++;
      }
   }
   MPI_Win_unlock(Proc->rank, Proc->Bound.Rma.win);
}

void DMEM_POS_AccumBoundAll(DMEM_ProcInfo *Proc,
                            SolveData *Out_data)
{
   int target;
   double start;
   for (int q = 0; q < Proc->Neighb.size; q++){
      target = Proc->Neighb.ranks[q];
      start = MPI_Wtime();
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
      Out_data->wtime_comm[0] += (MPI_Wtime() - start);
      for (int i = 0; i < Proc->Bound.Rma.send_count[q]; i++){
         Proc->Bound.Rma.send[q][i] = 0;
      }
   }
}

void DMEM_RMA_AccumBound(DMEM_ProcInfo *Proc,
                         int q)
{
   int target = Proc->Neighb.ranks[q];
   MPI_Accumulate(Proc->Bound.Rma.send[q],
                  Proc->Bound.Rma.send_count[q],
                  MPI_DOUBLE,
                  target,
                  Proc->Bound.Rma.send_disp[q],
                  Proc->Bound.Rma.send_count[q],
                  MPI_DOUBLE,
                  MPI_SUM,
                  Proc->Bound.Rma.win);
}

void DMEM_POS_AccumBound(DMEM_ProcInfo *Proc,
                         SolveData *Out_data,
                         int q)
{
   int target = Proc->Neighb.ranks[q];
   double start;
   start = MPI_Wtime();
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
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
   for (int i = 0; i < Proc->Bound.Rma.send_count[q]; i++){
      Proc->Bound.Rma.send[q][i] = 0;
   }
}

void DMEM_DiagSweep_GS(DMEM_MatrixInfo Mat,
                       DMEM_ProcInfo Proc,
                       SolveData *Out_data,
                       SolveVars *Vars)
{
   int ind, row, shift_row;
   double u_prev, elem, start;
   double relax = 0;

   start = MPI_Wtime();
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
   Out_data->flop_scaled[0] += (double)Mat.D.nnz/(double)Mat.n_glob;
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->relax_scaled[0] += (relax / (double)Mat.n_glob);
}

void DMEM_DiagSweep_DirectSolve(DMEM_MatrixInfo Mat,
                                DMEM_ProcInfo Proc,
                                SolveData *Out_data,
                                SolveVars *Vars)
{
   int ind, row, shift_row;
   double u_prev, elem, start;
   double relax = 0;
   start = MPI_Wtime();
   PARDISO(Mat.Pard.pt,
           &(Mat.Pard.maxfct),
           &(Mat.Pard.mnum),
           &(Mat.Pard.mtype),
           &(Mat.Pard.phase),
           &(Mat.Pard.csr.n),
           Mat.D.a,
           Mat.Pard.csr.ia,
           Mat.Pard.csr.ja,
           &(Mat.Pard.idum),
          // Mat.Pard.perm,
           &(Mat.Pard.nrhs),
           Mat.Pard.iparm,
           &(Mat.Pard.msglvl),
           Vars->r,
           Vars->du,
           &(Mat.Pard.error));
   for (int i = 0; i < Mat.n; i++){
      Vars->u[i] += Vars->du[i];
      Vars->r[i] = 0;
      relax++;
   }
   Out_data->relax[0]++;
   Out_data->flop_scaled[0] += 
      (double)Mat.Pard.iparm[17]/(double)Mat.n_glob;
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->relax_scaled[0] += (relax / (double)Mat.n_glob);
}

void DMEM_OffDiagSweep(DMEM_MatrixInfo Mat,
                       DMEM_ProcInfo *Proc,
                       SolveData *Out_data,
                       SolveVars *Vars)
{
   int ind, row, map, shift_row, k, neighb_rank, s;
   double elem, start;

   start = MPI_Wtime();
   for (int q = 0; q < Proc->Neighb.size; q++){
      neighb_rank = Proc->Neighb.ranks[q];
      for (int i = 0; i < Mat.n; i++){
         k = Mat.B[q].j_ptr[i];
         for (int j = 0; j < Mat.B[q].j_ptr[i+1]-Mat.B[q].j_ptr[i]; j++) {
            ind = k+j;
            elem = Mat.B[q].a[ind];
            row = Mat.B[q].i[ind];
            shift_row = row - Proc->Bound.row_min[q];
            map = Proc->Bound.send_points[q][shift_row];
            Proc->Bound.Rma.send[q][map] += (Vars->du[i] * elem);
         }
      }
   }
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->wtime_update_bound[0] += (MPI_Wtime() - start);
}

void DMEM_DiagSweep(DMEM_MatrixInfo Mat,
                    DMEM_ProcInfo *Proc,
                    SolveData *Out_data,
                    DMEM_SolveVars *Vars)
{
   if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
      DMEM_DiagSweep_DirectSolve(Mat, *Proc, Out_data, &(Vars->Loc));
   }
   else {
      DMEM_DiagSweep_GS(Mat, *Proc, Out_data, &(Vars->Loc));
   }
}

int DMEM_POS_CheckConverge(double norm, 
                           SolveData Out_data, 
                           SolveParams Params,
                           DMEM_ProcInfo *Proc)
{
   int one = 1;
   int read_recv;
   int target;
   if (((Out_data.sweep[0] >= Params.sweep_max) || (norm < Params.tol)) &&
       (!(Proc->Conv.my_flag))){
      for (int q = 0; q < Proc->Neighb.size; q++){
         target = Proc->Neighb.ranks[q];
         MPI_Win_lock(MPI_LOCK_SHARED, target, 0, Proc->Conv.win);
         MPI_Accumulate(&one, 1, MPI_INT, target, 0, 1,
                        MPI_INT, MPI_SUM, Proc->Conv.win);
         MPI_Win_unlock(target, Proc->Conv.win);
      }
      Proc->Conv.my_flag = 1;
   }
   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Conv.win);
   read_recv = Proc->Conv.recv[0];
   MPI_Win_unlock(Proc->rank, Proc->Conv.win);

   return read_recv;
}

int DMEM_SOS_CheckConverge(double norm,
                           DMEM_ProcInfo *Proc,
                           SolveData *Out_data,
                           SolveParams Params)
{
   int one = 1;
   int neighb_conv, conv = 0;
   int target;
   int ind;
 
   if (Out_data->sweep[0] >= Params.sweep_max) return 1;

   return 0;
}

int DMEM_POS_CheckConverge_array(double norm,
                                 SolveData Out_data,
                                 SolveParams Params,
                                 DMEM_ProcInfo *Proc)
{
   int one = 1;
   int read_recv;
   if (((Out_data.sweep[0] >= Params.sweep_max) || (norm < Params.tol)) &&
       (!(Proc->Conv.my_flag))){
      for (int p = 0; p < Proc->world_size; p++){
         MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p, 0, Proc->Conv.win);
         MPI_Put(&one, 1, MPI_INT, p, Proc->rank, 1,
                 MPI_INT, Proc->Conv.win);
         MPI_Win_unlock(p, Proc->Conv.win);
      }
      Proc->Conv.my_flag = 1;
   }

   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Conv.win);

   read_recv = SumInt(Proc->Conv.recv, Proc->world_size);

   MPI_Win_unlock(Proc->rank, Proc->Conv.win);
   return read_recv;
}

int DMEM_COLLECT_CheckConverge(int n,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data,
                               SolveParams Params)
{
   int flag = 0, sum_flags;
   double r_norm_glob, r_norm_squ, start;
  
   start = MPI_Wtime();
   r_norm_squ = SumSquaredDouble(Vars->Loc.r, n);
   MPI_Allreduce(&r_norm_squ,
                 &r_norm_glob,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 MPI_COMM_WORLD);
   r_norm_glob = sqrt(r_norm_glob);
   if ((r_norm_glob < Params.tol) ||
       (Out_data->sweep[0] >= Params.sweep_max)){
      flag = 1;
   }
   Out_data->wtime_conv[0] += (MPI_Wtime() - start);
   return flag;
}

void DMEM_SOS_AccumBoundAll(DMEM_ProcInfo *Proc, SolveData *Out_data)
{
   int target;
   double start;
   
   start = MPI_Wtime();
   MPI_Win_post(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
   MPI_Win_start(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
   for (int q = 0; q < Proc->Neighb.size; q++){
      target = Proc->Neighb.ranks[q];
      MPI_Accumulate(Proc->Bound.Rma.send[q],
                     Proc->Bound.Rma.send_count[q],
                     MPI_DOUBLE,
                     target,
                     Proc->Bound.Rma.send_disp[q],
                     Proc->Bound.Rma.send_count[q],
                     MPI_DOUBLE,
                     MPI_SUM,
                     Proc->Bound.Rma.win);
   }
   MPI_Win_complete(Proc->Bound.Rma.win);
   MPI_Win_wait(Proc->Bound.Rma.win);
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
   Out_data->wtime_comm_sweep[0] += (MPI_Wtime() - start);
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int i = 0; i < Proc->Bound.Rma.send_count[q]; i++){
         Proc->Bound.Rma.send[q][i] = 0;
      }
   }
}

void DMEM_SOS_Recv(DMEM_ProcInfo *Proc,
                   SolveData *Out_data)
{
   double start;

   start = MPI_Wtime();
   MPI_Win_post(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
   MPI_Win_start(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
   MPI_Win_complete(Proc->Bound.Rma.win);
   MPI_Win_wait(Proc->Bound.Rma.win);
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
}
