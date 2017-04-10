#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"


void DMEM_RMA_MCGS_RecvBound(DMEM_MulticolorInfo *MC_info,
                             DMEM_ProcInfo Proc,
                             double *r,
                             int c)
{
   int map, recv_disp;
   double elem;
   int disp_loc = MC_info->disp_loc[c];

   int k = 0, j = 0;
   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, Proc.rank, 0, MC_info->Rma[c].win);
   for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
      for (int q = 0; q < MC_info->NeighbRecv[c][cc].size; q++){
         for (int i = 0; i < MC_info->Bound[c][cc].num_recv_points[q]; i++){
            map = MC_info->Bound[c][cc].recv_points[q][i];
            elem = MC_info->Rma[c].recv[k];
            r[map] -= elem;
            MC_info->Rma[c].recv[k] = 0;
            k++;
         }
      }
   }
   MPI_Win_unlock(Proc.rank, MC_info->Rma[c].win);
}

void DMEM_MCGS_DiagSweep_GS(DMEM_MulticolorInfo MC_info,
                            DMEM_ProcInfo Proc,
                            SolveData *Out_data,
                            SolveVars *Vars,
                            int c)
{
   int ind, row, shift_row;
   double u_prev, elem, start;
   double relax = 0;
   int r_disp_loc = MC_info.disp_loc[c];
   int disp_glob =  MC_info.P_glob.disp[c];
   int disp_loc = MC_info.P_loc[c].disp[Proc.rank];

   start = MPI_Wtime();
   for (int i = 0; i < MC_info.D[c].n; i++){
      u_prev = Vars->u[r_disp_loc+i];
      Vars->u[r_disp_loc+i] += 
         (Vars->r[r_disp_loc+i] / MC_info.D[c].diag[i]);
      Vars->du[r_disp_loc+i] = Vars->u[r_disp_loc+i] - u_prev;
      for (int j = 0; j < 
              MC_info.D[c].j_ptr[i+1]-MC_info.D[c].j_ptr[i]; j++) {
         ind = MC_info.D[c].j_ptr[i]+j;
         row = MC_info.D[c].i[ind];
         elem = MC_info.D[c].a[ind];
         shift_row = row - (disp_glob + disp_loc);
         Vars->r[r_disp_loc+shift_row] -=
            (Vars->du[r_disp_loc+i] * elem);
      }
      relax++;
      Out_data->relax[0]++;
   }
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
   Out_data->relax_scaled[0] += (relax / (double)MC_info.n_glob);
}

void DMEM_MCGS_OffDiagSweep_GS(DMEM_MulticolorInfo MC_info,
                               DMEM_ProcInfo Proc,
                               SolveData *Out_data,
                               SolveVars *Vars,
                               int c)
{
   int ind, row, map, shift_row, k, neighb_rank, s;
   double elem, start;
   int r_disp_loc = MC_info.disp_loc[c];

   start = MPI_Wtime();
  
   for (int cc = 0; cc < MC_info.P_glob.nparts; cc++){ 
      for (int q = 0; q < MC_info.NeighbSend[c][cc].size; q++){
         neighb_rank = MC_info.NeighbSend[c][cc].ranks[q];
         for (int i = 0; i < MC_info.D[c].n; i++){
            k = MC_info.B[c][cc][q].j_ptr[i];
            for (int j = 0; j < MC_info.B[c][cc][q].j_ptr[i+1] - 
                                MC_info.B[c][cc][q].j_ptr[i]; j++) {
               ind = k+j;
               elem = MC_info.B[c][cc][q].a[ind];
               row = MC_info.B[c][cc][q].i[ind];
               shift_row = row - MC_info.Bound[c][cc].row_min[q];
               map = MC_info.Bound[c][cc].send_points[q][shift_row];
               MC_info.Bound[c][cc].Rma.send[q][map] += 
                  (Vars->du[r_disp_loc+i] * elem);
            }
         }
      }
   }
   Out_data->wtime_comp[0] += (MPI_Wtime() - start);
}

void DMEM_SOS_MCGS_AccumBoundAll(DMEM_MulticolorInfo *MC_info, 
                                 DMEM_ProcInfo Proc,
                                 SolveData *Out_data,
                                 int c)
{
   int target;
   double start;

   start = MPI_Wtime();
   MPI_Win_post(MC_info->mpi_group, 0, MC_info->Rma[c].win);
   MPI_Win_start(MC_info->mpi_group, 0, MC_info->Rma[c].win);
   for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
      for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
         target = MC_info->NeighbSend[c][cc].ranks[q];
         MPI_Accumulate(MC_info->Bound[c][cc].Rma.send[q],
                        MC_info->Bound[c][cc].Rma.send_count[q],
                        MPI_DOUBLE,
                        target,
                        MC_info->Bound[c][cc].Rma.send_disp[q],
                        MC_info->Bound[c][cc].Rma.send_count[q],
                        MPI_DOUBLE,
                        MPI_SUM,
                        MC_info->Rma[c].win);
      }
   }
   MPI_Win_complete(MC_info->Rma[c].win);
   MPI_Win_wait(MC_info->Rma[c].win);
   Out_data->wtime_comm[0] += (MPI_Wtime() - start);
   for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
      for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
         for (int i = 0; i < MC_info->Bound[c][cc].Rma.send_count[q]; i++){
            MC_info->Bound[c][cc].Rma.send[q][i] = 0;
         }
      }
   }
}

void DMEM_MC_ZeroData(DMEM_MulticolorInfo *MC_info)
{
   for (int c = 0; c < MC_info->P_glob.nparts; c++){
      for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
         for (int q = 0; q < MC_info->NeighbSend[c][cc].size; q++){
            for (int i = 0; i < 
                    MC_info->Bound[c][cc].Rma.send_count[q]; i++){
               MC_info->Bound[c][cc].Rma.send[q][i] = 0;
            }
         }
      }
      for (int i = 0; i < MC_info->Rma[c].recv_count; i++){
         MC_info->Rma[c].recv[i] = 0;
      }
   }
}
