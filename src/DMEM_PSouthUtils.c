#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "SEQ_Solve.h"

void DMEM_RMA_PSOUTH_RecvResNorm(DMEM_ProcInfo *Proc)
{
   double elem, start;
   int k = 0;

  // MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Res.Rma.win);
   for (int q = 0; q < Proc->Neighb.size; q++){
      if (!AlmostEqual2sComplement_DBL(Proc->Res.Rma.recv[q],
                                       Proc->Res.Rma.recv_prev[q], 1)){
         Proc->Res.norm[q+1] = Proc->Res.Rma.recv[q];
         Proc->Res.Rma.recv_prev[q] =  Proc->Res.Rma.recv[q];
      }
   }
  // MPI_Win_unlock(Proc->rank, Proc->Res.Rma.win);
}

void DMEM_RMA_PSOUTH_RecvBound(DMEM_MatrixInfo Mat,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars)
{
   int map;
   double elem, diff, recv_norm;
   int k = 0;

  // MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc->rank, 0, Proc->Bound.Rma.win);
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int j = 0; j < Proc->Bound.num_recv_points[q]; j++){
         map = Proc->Bound.recv_points[q][j];
         Vars->Loc.r[map] -= Proc->Bound.Rma.recv[k];
         Proc->Bound.Rma.recv[k] = 0;
         k++;
      }
     // if (!AlmostEqual2sComplement_DBL(Proc->Bound.Rma.recv[k],
     //                                  Proc->Bound.Rma.recv_prev[k], 1)){
         Proc->Res.norm[q+1] =  Proc->Bound.Rma.recv[k];
         Proc->Bound.Rma.recv_prev[k] =  Proc->Bound.Rma.recv[k];
     // }
      k++;
   }
  // MPI_Win_unlock(Proc->rank, Proc->Bound.Rma.win);
}

void DMEM_SOS_PSOUTH_PutResNormAll(DMEM_ProcInfo *Proc,
                                   SolveData *Out_data)
{
   int target;
   double start;

   start = MPI_Wtime();
   MPI_Win_post(Proc->Neighb.mpi_group, 0, Proc->Res.Rma.win);
   MPI_Win_start(Proc->Neighb.mpi_group, 0, Proc->Res.Rma.win);
   for (int q = 0; q < Proc->Neighb.size; q++){
      target = Proc->Neighb.ranks[q];
      MPI_Put(&(Proc->Res.norm[0]),
              1,
              MPI_DOUBLE,
              target,
              Proc->Res.Rma.send_disp[q],
              1,
              MPI_DOUBLE,
              Proc->Res.Rma.win);
   }
   MPI_Win_complete(Proc->Res.Rma.win);
   MPI_Win_wait(Proc->Res.Rma.win);
   Out_data->comm[0] += Proc->Neighb.size;
   Out_data->comm_scaled[0] +=
         ((double)Proc->Neighb.size/(double)Proc->world_size);
   Out_data->comm_res[0] += Proc->Neighb.size;
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
   Out_data->wtime_comm_res[0] += (MPI_Wtime() - start);
}
