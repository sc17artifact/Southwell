#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "DMEM_SolveUtils.h"
#include "SEQ_Solve.h"

void DMEM_DistrSouthwellAsync_Pos(DMEM_MatrixInfo Mat,
                                  SolveParams Params,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveData *Out_data)
{
   int k, row, col, shift_row, map, conv_flag, update_flag;
   double elem, r_norm, start;
   if (Proc->world_size == 1) {
      SEQ_JacobiRes(Mat.A, &(Vars->Loc), Params, Out_data);
   }
   else {
      start = MPI_Wtime();
      while(1){
         update_flag = IsZeroMaxAbsDouble(Proc->Res.norm1, 
                                          Proc->Neighb.size+1);
         if (update_flag){
            RMA_Recv(Mat, Proc, Vars->Loc.r);
            DiagSweep_GS(Mat, *Proc, Out_data, &(Vars->Loc));
            for (int q = 0; q < Proc->Neighb.size; q++){
               for (int i = 0; i < Proc->Bound.unique_size[q]; i++){
                  Proc->Bound.Rma.send_prev[q][i] = 
                     Proc->Bound.Rma.send[q][i];
               }
            }
            OffDiagSweep_GS(Mat, Proc, Out_data, &(Vars->Loc));
           // DMEM_UpdateResEstim_Norm1(Mat, Proc);
            POS_SendAll(Proc);
            Proc->Res.norm1[0] = r_norm = Norm1(Vars->Loc.r, Mat.n);
         }
         Out_data->sweep[0]++;
         conv_flag = POS_CheckConverge(r_norm, *Out_data, Params, Proc);
         if (conv_flag == Proc->world_size){
            break;
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      Out_data->wtime_tot[0] = MPI_Wtime() - start;
   }
}

void DMEM_JacobiAsync_Pos(DMEM_MatrixInfo Mat,
                          SolveParams Params,
                          DMEM_ProcInfo *Proc,
                          DMEM_SolveVars *Vars,
                          SolveData *Out_data)
{ 
   int k, row, col, shift_row, map, conv_flag;
   double elem, r_norm, start;
   if (Proc->world_size == 1) {
      SEQ_JacobiRes(Mat.A, &(Vars->Loc), Params, Out_data);
   }
   else {
      start = MPI_Wtime();
      while(1){
         RMA_Recv(Mat, Proc, Vars->Loc.r);
         DiagSweep_GS(Mat, *Proc, Out_data, &(Vars->Loc));
         OffDiagSweep_GS(Mat, Proc, Out_data, &(Vars->Loc));
         POS_SendAll(Proc);
         r_norm = Norm1(Vars->Loc.r, Mat.n);
         Out_data->sweep[0]++;
         conv_flag = POS_CheckConverge(r_norm, *Out_data, Params, Proc);
         if (conv_flag == Proc->world_size){
            break;
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      Out_data->wtime_tot[0] = MPI_Wtime() - start;
   }
}
