#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "DMEM_SolveUtils.h"
#include "DMEM_PSouthUtils.h"
#include "DMEM_DSouthUtils.h"
#include "DMEM_MCGSutils.h"
#include "SEQ_Solve.h"

int ds_fail = 0, ds_pass = 0;

void DMEM_POS_AsyncDistrSouthwell(DMEM_MatrixInfo Mat,
                                  SolveParams Params,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveData *Out_data)
{
//   int row, col, shift_row, map;
//   int neighb_conv_flag, update_flag;
//   int k;
//   double elem, r_norm, start;
//   Proc->Res.explicit_flag = 1;
//   Proc->Conv.my_flag = 0;
//   VSLStreamStatePtr mkl_rand_stream;
//   vslNewStream(&mkl_rand_stream, VSL_BRNG_NONDETERM, time(NULL)*Proc->rank);
//
//   if ((Params.sweep_max == 0) ||
//       (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
//      return;
//   }
//  // Proc->Bound.res_norm = sqrt(DMEM_BoundSumSqu(Proc, *Vars));
//   Proc->Res.norm[0] = 
//      Proc->Bound.res_norm = 
//         Norm2(Vars->Loc.r, Mat.n);
//   for (int q = 0; q < Proc->Neighb.size; q++){
//      DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
//   }
//   DMEM_SOS_AccumBoundAll(Proc, Out_data);
//   MPI_Barrier(MPI_COMM_WORLD);
//   DMEM_RMA_DSOUTH_RecvAllResBound(Mat, Proc, Vars, Params);
//
//   for (int q = 0; q < Proc->Neighb.size; q++){
//      DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
//   }
//   DMEM_SOS_AccumBoundAll(Proc, Out_data);
//   MPI_Barrier(MPI_COMM_WORLD);
//   DMEM_RMA_DSOUTH_RecvAllResBound(Mat, Proc, Vars, Params);
//
//   Vars->Loc.ds_delay_count = Params.ds_delay; 
//   start = MPI_Wtime();
//   while(1){
//     // Proc->Bound.res_norm = sqrt(DMEM_BoundSumSqu(Proc, *Vars));
//      update_flag = IsZeroMaxAbsDouble(Proc->Res.norm,
//                                       Proc->Neighb.size+1);
//      if (update_flag){
//         Vars->Loc.ds_delay_count = 0;
//         DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
//         DMEM_DSOUTH_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
//         if (ds_res_estim_flag){
//            DMEM_DSOUTH_UpdateResEstim_Norm2(Mat, Proc, Out_data);
//         }
//         if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
//            Proc->Res.norm[0] = 
//               Proc->Res.norm_prev[0] = 
//                  Proc->Bound.res_norm = 0;
//         }
//         else {
//            Proc->Res.norm[0] = 
//               Proc->Res.norm_prev[0] = 
//                  Proc->Bound.res_norm = 
//                     Norm2(Vars->Loc.r, Mat.n);
//         }
//        // Proc->Bound.res_norm = sqrt(DMEM_BoundSumSqu(Proc, *Vars));
//         for (int q = 0; q < Proc->Neighb.size; q++){
//            DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
//         }
//         DMEM_POS_AccumBoundAll(Proc, Out_data);
//         Proc->Res.explicit_flag = 0;
//         for (int q = 0; q < Proc->Neighb.size; q++){
//            Proc->Res.neighb_explicit_flag[q] = 0;
//         }
//         Out_data->comm[0] += Proc->Neighb.size;
//         Out_data->comm_scaled[0] += 
//            ((double)Proc->Neighb.size/(double)Proc->world_size);
//         Out_data->comm_sweep[0] += Proc->Neighb.size;
//      }
//      DMEM_RMA_DSOUTH_RecvAllBound(Mat, Proc, Vars, Params);
//      Proc->Res.norm[0] =
//         Proc->Bound.res_norm =
//            Norm2(Vars->Loc.r, Mat.n);
//      Out_data->sweep[0]++;
//      neighb_conv_flag = DMEM_POS_CheckConverge(Proc->Res.norm[0],
//                                                *Out_data,
//                                                Params,
//                                                Proc);
//      if ((neighb_conv_flag == Proc->Neighb.size) &&
//          (Proc->Conv.my_flag == 1)){
//         break;
//      }
//      if (!AlmostEqual2sComplement_DBL(Proc->Res.norm[0], 
//                                       Proc->Res.norm_prev[0], 1)){
//         Vars->Loc.ds_delay_count = 0;
//         Proc->Res.norm_prev[0] = Proc->Res.norm[0];
//         Proc->Res.explicit_flag = 1;
//         for (int q = 0; q < Proc->Neighb.size; q++){
//            Proc->Res.neighb_explicit_flag[q] = 1;
//         }
//      }
//      if (Vars->Loc.ds_delay_count == Params.ds_delay){
//         DMEM_DSOUTH_ExplicitResUpdate(Mat,
//                                       Proc,
//                                       Out_data,
//                                       Vars,
//                                       Params,
//                                       mkl_rand_stream);
//         Vars->Loc.ds_delay_count = 0;
//      }
//      else {
//         Vars->Loc.ds_delay_count++;
//      }
//      DMEM_RMA_DSOUTH_RecvAllBound(Mat, Proc, Vars, Params);
//   }
//   Out_data->wtime_tot[0] = MPI_Wtime() - start; 
//  // for (int p = 0; p < Proc->world_size; p++){
//  //    if (p == Proc->rank){
//  //       printf("(%d, %e): ", p, Proc->Res.norm[0]);
//  //       for (int q = 0; q < Proc->Neighb.size; q++){
//  //          printf(" (%d, %e, %e) ", Proc->Neighb.ranks[q], Proc->Res.neighb_norm_estim[q], Proc->Res.norm[q+1]);
//  //       }
//  //       printf("\n");
//  //    }
//  //    MPI_Barrier(MPI_COMM_WORLD);
//  // }
//  // for (int p = 0; p < Proc->world_size; p++){
//  //    if (p == Proc->rank){
//  //       printf("(%d, %e): ", p, Proc->Res.norm[0]);
//  //       for (int q = 0; q < Proc->Neighb.size; q++){
//  //          printf(" (%d, %e) ", Proc->Neighb.ranks[q], Proc->Res.neighb_norm_estim[q]);
//  //       }
//  //       printf("\n");
//  //    }
//  //    MPI_Barrier(MPI_COMM_WORLD);
//  // }
}


void DMEM_SOS_SyncDistrSouthwell(DMEM_MatrixInfo Mat,
                                 SolveParams Params,
                                 DMEM_ProcInfo *Proc,
                                 DMEM_SolveVars *Vars,
                                 SolveData *Out_data)
{
   int row, col, shift_row, map;
   int conv_flag, update_flag;
   int k;
   double start_iter, elem, r_norm, start, prob;
   VSLStreamStatePtr mkl_rand_stream;
   vslNewStream(&mkl_rand_stream, VSL_BRNG_NONDETERM, time(NULL)*Proc->rank);

   if ((Params.sweep_max == 0) ||
       (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
      return;
   }
   Proc->Res.norm[0] =
      Proc->Bound.res_norm =
         Norm2(Vars->Loc.r, Mat.n);
   for (int q = 0; q < Proc->Neighb.size; q++){
      DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
   }
   DMEM_SOS_AccumBoundAll(Proc, Out_data);
   MPI_Barrier(MPI_COMM_WORLD);
   DMEM_RMA_DSOUTH_RecvAllResBound(Mat, Proc, Vars, Params);

   for (int q = 0; q < Proc->Neighb.size; q++){
      DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
   }
   DMEM_SOS_AccumBoundAll(Proc, Out_data);
   MPI_Barrier(MPI_COMM_WORLD);
   DMEM_RMA_DSOUTH_RecvAllResBound(Mat, Proc, Vars, Params);
   MPI_Barrier(MPI_COMM_WORLD);

   Out_data->work_iter[0] = 0;
   Vars->Loc.ds_delay_count = Params.ds_delay;
  // if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
  //    Out_data->wtime_tot[0] = Mat.Pard.wtime_setup;
  //    Out_data->wtime_comp[0] = Mat.Pard.wtime_setup;
  // }
  // else {
      Out_data->wtime_tot[0] = 0;
  // }
   start = MPI_Wtime();
   while(1){
       update_flag = IsZeroMaxAbsDouble(Proc->Res.norm,
                                        Proc->Neighb.size+1);
      
      Out_data->relax_mask[0] = 0;
      if (update_flag){
         Out_data->relax_mask[0] = 1;
         Vars->Loc.ds_delay_count = 0;
         DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
         DMEM_DSOUTH_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
         if (ds_res_estim_flag){
            DMEM_DSOUTH_UpdateResEstim_Norm2(Mat, Proc, Out_data);
         }
         if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
            Proc->Res.norm[0] =
               Proc->Res.norm_prev[0] =
                  Proc->Bound.res_norm = 0;
         }
         else {
            Proc->Res.norm[0] =
               Proc->Res.norm_prev[0] =
                  Proc->Bound.res_norm =
                     Norm2(Vars->Loc.r, Mat.n);
         }
         for (int q = 0; q < Proc->Neighb.size; q++){
            DMEM_RMA_DSOUTH_GatherBoundRes(Mat, Proc, *Vars, Out_data, q);
         }
         DMEM_SOS_AccumBoundAll(Proc, Out_data);
         Proc->Res.explicit_flag = 0;
         for (int q = 0; q < Proc->Neighb.size; q++){
            Proc->Res.neighb_explicit_flag[q] = 0;
         }
         Out_data->comm[0] += Proc->Neighb.size;
         Out_data->comm_scaled[0] +=
            ((double)Proc->Neighb.size/(double)Proc->world_size);
         Out_data->comm_sweep[0] += Proc->Neighb.size;
      }
      else {
         DMEM_SOS_Recv(Proc, Out_data);
      }
      DMEM_RMA_DSOUTH_RecvAllBound(Mat, Proc, Vars, Params);
      Proc->Res.norm[0] =
         Proc->Bound.res_norm =
            Norm2(Vars->Loc.r, Mat.n);

      Out_data->sweep[0]++;
      conv_flag = DMEM_SOS_CheckConverge(Proc->Res.norm[0],
                                         Proc,
                                         Out_data,
                                         Params);
      if (conv_flag) break;
      if (!AlmostEqual2sComplement_DBL(Proc->Res.norm[0], 
                                       Proc->Res.norm_prev[0], 1)){
         Vars->Loc.ds_delay_count = 0;
         Proc->Res.norm_prev[0] = Proc->Res.norm[0];
         Proc->Res.explicit_flag = 1;
         for (int q = 0; q < Proc->Neighb.size; q++){
            Proc->Res.neighb_explicit_flag[q] = 1;
         }
      }
      start_iter = MPI_Wtime();
      MPI_Win_post(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
      MPI_Win_start(Proc->Neighb.mpi_group, 0, Proc->Bound.Rma.win);
      if (Vars->Loc.ds_delay_count == Params.ds_delay){
         DMEM_DSOUTH_ExplicitResUpdate(Mat,
                                       Proc,
                                       Out_data,
                                       Vars,
                                       Params,
                                       mkl_rand_stream);
         Vars->Loc.ds_delay_count = 0;
      }
      else {
         Vars->Loc.ds_delay_count++;
      }
      MPI_Win_complete(Proc->Bound.Rma.win);
      MPI_Win_wait(Proc->Bound.Rma.win);
      Out_data->wtime_comm_res[0] += (MPI_Wtime() - start_iter);
      Out_data->wtime_comm[0] += (MPI_Wtime() - start_iter);
      DMEM_RMA_DSOUTH_RecvAllBound(Mat, Proc, Vars, Params);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   Out_data->wtime_tot[0] += (MPI_Wtime() - start);
  // printf("%d CONVERGED!\n", Proc->rank);
  // for (int p = 0; p < Proc->world_size; p++){
  //    if (p == Proc->rank){
  //       printf("%d %d\n", p, Out_data->relax_mask[0]);
  //    }
  //    MPI_Barrier(MPI_COMM_WORLD);
  // }
  // for (int p = 0; p < Proc->world_size; p++){
  //    if (p == Proc->rank){
  //       printf("(%d, %e): ", p, Proc->Res.norm[0]);
  //       for (int q = 0; q < Proc->Neighb.size; q++){
  //          printf(" (%d, %e) ", Proc->Neighb.ranks[q], Proc->Res.norm[q+1]);
  //       }
  //       printf("\n");
  //    }
  //    MPI_Barrier(MPI_COMM_WORLD);
  // }
  // for (int p = 0; p < Proc->world_size; p++){
  //    if (p == Proc->rank){
  //       printf("(%d, %e): ", p, Proc->Res.norm[0]);
  //       for (int q = 0; q < Proc->Neighb.size; q++){
  //          printf(" (%d, %e) ", 
  //                 Proc->Neighb.ranks[q], 
  //                 Proc->Res.neighb_norm_estim[q]);
  //       }
  //       printf("\n");
  //    }
  //    MPI_Barrier(MPI_COMM_WORLD);
  // }
}


void DMEM_POS_AsyncParSouthwell(DMEM_MatrixInfo Mat,
                                SolveParams Params,
                                DMEM_ProcInfo *Proc,
                                DMEM_SolveVars *Vars,
                                SolveData *Out_data)
{
//   int row, col, shift_row, map, end;
//   int neighb_conv_flag, update_flag;
//   int k;
//   double elem, r_norm, start;
//
//   if ((Params.sweep_max == 0) ||
//       (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
//      return;
//   }
//   Proc->Res.norm[0]  = Norm2(Vars->Loc.r, Mat.n);
//   for (int q = 0; q < Proc->Neighb.size; q++){
//      end = Proc->Bound.num_send_points[q];
//      Proc->Bound.Rma.send[q][end] = Proc->Res.norm[0];
//   }
//   Proc->Res.norm_prev[0] = Proc->Res.norm[0];
//   DMEM_POS_AccumBoundAll(Proc, Out_data);
//   MPI_Barrier(MPI_COMM_WORLD);
//   start = MPI_Wtime();
//   while(1){
//      DMEM_RMA_PSOUTH_RecvBound(Mat, Proc, Vars, Params);
//      Proc->Res.norm[0]  = Norm2(Vars->Loc.r, Mat.n);
//      update_flag = IsZeroMaxAbsDouble(Proc->Res.norm,
//                                       Proc->Neighb.size+1);
//      if (update_flag){
//         DMEM_RMA_PSOUTH_RecvBound(Mat, Proc, Vars, Params);
//         DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
//         DMEM_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
//         DMEM_RMA_PSOUTH_RecvBound(Mat, Proc, Vars);
//         if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
//            Proc->Res.norm[0] =
//               Proc->Res.norm_prev[0] =
//                  Proc->Bound.res_norm = 0;
//         }
//         else {
//            Proc->Res.norm[0] =
//               Proc->Res.norm_prev[0] =
//                  Proc->Bound.res_norm =
//                     Norm2(Vars->Loc.r, Mat.n);
//         }
//         for (int q = 0; q < Proc->Neighb.size; q++){
//            end = Proc->Bound.num_send_points[q];
//            Proc->Bound.Rma.send[q][end] = 
//              (Proc->Res.norm[0] - Proc->Res.norm_prev[0]);
//         }
//         DMEM_POS_AccumBoundAll(Proc, Out_data);
//         Out_data->comm[0] += Proc->Neighb.size;
//         Out_data->comm_scaled[0] +=
//            ((double)Proc->Neighb.size/(double)Proc->world_size);
//      }
//      Out_data->sweep[0]++;
//      neighb_conv_flag = DMEM_POS_CheckConverge(Proc->Res.norm[0],
//                                                *Out_data,
//                                                Params,
//                                                Proc);
//      if ((neighb_conv_flag == Proc->Neighb.size) &&
//          (Proc->Conv.my_flag == 1)){
//         break;
//      }
//   }
//   Out_data->wtime_tot[0] = MPI_Wtime() - start;
//
// // for (int p = 0; p < Proc->world_size; p++){
// //    if (p == Proc->rank){
// //       printf("(%d, %e): ", p, Proc->Res.norm[0]);
// //       for (int q = 0; q < Proc->Neighb.size; q++){
// //          printf(" (%d, %e) ", Proc->Neighb.ranks[q], Proc->Res.norm[q+1]);
// //       }
// //       printf("\n");
// //    }
// //    MPI_Barrier(MPI_COMM_WORLD);
// // }

}

void DMEM_SOS_SyncParSouthwell(DMEM_MatrixInfo Mat,
                               SolveParams Params,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data)
{
   int k, row, col, shift_row, map; 
   int update_flag = 0, conv_flag, end;
   double elem, r_norm, start;
   if ((Params.sweep_max == 0) ||
       (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
      return;
   }
   Out_data->work_iter[0] = 0;
  // if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
  //    Out_data->wtime_tot[0] = Mat.Pard.wtime_setup;
  //    Out_data->wtime_comp[0] = Mat.Pard.wtime_setup;
  // }
  // else {
      Out_data->wtime_tot[0] = 0;
  // }
   start = MPI_Wtime();
   while(1){
      Proc->Res.norm_prev[0] =  Proc->Res.norm[0];
      Proc->Res.norm[0] = Norm2(Vars->Loc.r, Mat.n);
      if (!AlmostEqual2sComplement_DBL(Proc->Res.norm[0],
                                       Proc->Res.norm_prev[0], 1)){
         DMEM_SOS_PSOUTH_PutResNormAll(Proc, Out_data);
      }
      else {
         DMEM_SOS_Recv(Proc, Out_data);
      }
      DMEM_RMA_PSOUTH_RecvResNorm(Proc);
      update_flag = IsZeroMaxAbsDouble(Proc->Res.norm,
                                       Proc->Neighb.size+1);
      Out_data->relax_mask[0] = 0;
      if (update_flag){ 
         Out_data->relax_mask[0] = 1;
         DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
         DMEM_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
         if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
            Proc->Res.norm[0] = 
               Proc->Res.norm_prev[0] = 0;
         }
         else {
            Proc->Res.norm[0] = 
               Proc->Res.norm_prev[0] = 
                  Norm2(Vars->Loc.r, Mat.n);
         }
         for (int q = 0; q < Proc->Neighb.size; q++){
            end = Proc->Bound.num_send_points[q];
            Proc->Bound.Rma.send[q][end] =
              (Proc->Res.norm[0] - Proc->Res.norm_send_prev[0]);
         }
         Proc->Res.norm_send_prev[0] = Proc->Res.norm[0];
         DMEM_SOS_AccumBoundAll(Proc, Out_data);
         Out_data->comm[0] += Proc->Neighb.size;
         Out_data->comm_scaled[0] +=
            ((double)Proc->Neighb.size/(double)Proc->world_size);
         Out_data->comm_sweep[0] += Proc->Neighb.size;
      }
      else {
         DMEM_SOS_Recv(Proc, Out_data);
      }
      DMEM_RMA_PSOUTH_RecvBound(Mat, Proc, Vars);
      Out_data->sweep[0]++;
      conv_flag = DMEM_SOS_CheckConverge(Proc->Res.norm[0],
                                         Proc,
                                         Out_data,
                                         Params);
      if (conv_flag) break;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   Out_data->wtime_tot[0] += (MPI_Wtime() - start);
  // Proc->Res.norm_prev[0] =  Proc->Res.norm[0];
  // Proc->Res.norm[0] = Norm2(Vars->Loc.r, Mat.n);
  // if (!AlmostEqual2sComplement_DBL(Proc->Res.norm[0],
  //                                  Proc->Res.norm_prev[0], 1)){
  //    DMEM_SOS_PSOUTH_PutResNormAll(Proc, Out_data);
  // }
  // else {
  //    DMEM_SOS_Recv(Proc, Out_data);
  // }
  // DMEM_RMA_PSOUTH_RecvResNorm(Proc);
  // update_flag = IsZeroMaxAbsDouble(Proc->Res.norm,
  //                                  Proc->Neighb.size+1);
  // for (int p = 0; p < Proc->world_size; p++){
  //    if (p == Proc->rank){
  //       printf("(%d, %e, %d):\n\t", p, Proc->Res.norm[0], update_flag);
  //       for (int q = 0; q < Proc->Neighb.size; q++){
  //          printf("(%d, %e) ", Proc->Neighb.ranks[q], Proc->Res.norm[q+1]);
  //       }
  //       printf("\n");
  //    }
  //    MPI_Barrier(MPI_COMM_WORLD);
  // }
}

void DMEM_POS_AsyncBlockJacobi(DMEM_MatrixInfo Mat,
                               SolveParams Params,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data)
{
  // int k, row, col, shift_row, map, neighb_conv_flag;
  // double elem, r_norm, start;
  // Proc->Conv.my_flag = 0;
  // if ((Params.sweep_max == 0) ||
  //     (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
  //    return;
  // }
  // Out_data->work_iter[0] = 1;
  // start = MPI_Wtime();
  // while(1){
  //    DMEM_RMA_RecvBound(Mat, Proc, Vars->Loc.r);
  //    DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
  //    DMEM_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
  //    DMEM_POS_AccumBoundAll(Proc, Out_data);
  //    Proc->Res.norm[0]  = Norm2(Vars->Loc.r, Mat.n);
  //    Out_data->comm[0] += Proc->Neighb.size;
  //    Out_data->comm_scaled[0] +=
  //          ((double)Proc->Neighb.size/(double)Proc->world_size);
  //    Out_data->sweep[0]++;
  //    neighb_conv_flag = DMEM_POS_CheckConverge(Proc->Res.norm[0],
  //                                              *Out_data,
  //                                              Params,
  //                                              Proc);
  //    if ((neighb_conv_flag == Proc->Neighb.size) &&
  //        (Proc->Conv.my_flag == 1)){
  //       break;
  //    }
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // Out_data->wtime_tot[0] = MPI_Wtime() - start;
}

void DMEM_SOS_SyncBlockJacobi(DMEM_MatrixInfo Mat,
                              SolveParams Params,
                              DMEM_ProcInfo *Proc,
                              DMEM_SolveVars *Vars,
                              SolveData *Out_data)
{
   int k, row, col, shift_row, map, conv_flag;
   double elem, start;
   if ((Params.sweep_max == 0) ||
       (Norm2(Vars->Glob.r, Mat.n_glob) < Params.tol)){
      return;
   }
   Out_data->work_iter[0] = 1;
   Out_data->relax_mask[0] = 1;
  // if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
  //    Out_data->wtime_tot[0] = Mat.Pard.wtime_setup;
  //    Out_data->wtime_comp[0] = Mat.Pard.wtime_setup;
  // }
  // else {
      Out_data->wtime_tot[0] = 0;
  // }
   start = MPI_Wtime();
   while(1){
      DMEM_DiagSweep(Mat, Proc, Out_data, Vars);
      DMEM_OffDiagSweep(Mat, Proc, Out_data, &(Vars->Loc));
      DMEM_SOS_AccumBoundAll(Proc, Out_data);
      MPI_Barrier(MPI_COMM_WORLD);
      DMEM_RMA_RecvBound(Mat, Proc, Vars->Loc.r);
      Out_data->comm[0] += Proc->Neighb.size;
      Out_data->comm_scaled[0] +=
            ((double)Proc->Neighb.size/(double)Proc->world_size);
      Out_data->sweep[0]++;
      conv_flag = DMEM_SOS_CheckConverge(Proc->Res.norm[0],
                                         Proc,
                                         Out_data,
                                         Params);
      if (conv_flag) break;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   Out_data->wtime_tot[0] += (MPI_Wtime() - start); 
}

void DMEM_SOS_MulticolorGS(DMEM_MulticolorInfo *MC_info,
                           DMEM_ProcInfo Proc,
                           SolveParams Params,
                           DMEM_SolveVars *Vars,
                           SolveData *Out_data)
{
   int k, row, col, shift_row, map, disp_loc, disp_glob;
   double elem, start;
   int conv_flag = 0;
   if ((Params.sweep_max == 0) ||
       (Norm2(Vars->Glob.r, MC_info->n_glob) < Params.tol)){
      return;
   }
   Out_data->work_iter[0] = 1;
   Out_data->relax_mask[0] = 1;
   start = MPI_Wtime();
   while(1){
      for (int c = 0; c < MC_info->P_glob.nparts; c++){
         DMEM_MCGS_DiagSweep_GS(*MC_info,
                                Proc,
                                Out_data,
                                &(Vars->Loc),
                                c);
         DMEM_MCGS_OffDiagSweep_GS(*MC_info,
                                   Proc,
                                   Out_data,
                                   &(Vars->Loc),
                                   c);
         DMEM_SOS_MCGS_AccumBoundAll(MC_info, Proc, Out_data, c);
         for (int cc = 0; cc < MC_info->P_glob.nparts; cc++){
            Out_data->comm[0] += MC_info->NeighbSend[c][cc].size;
            Out_data->comm_scaled[0] +=
               ((double)MC_info->NeighbSend[c][cc].size/
                               (double)Proc.world_size);
         }
         DMEM_RMA_MCGS_RecvBound(MC_info,
                                 Proc,
                                 Vars->Loc.r,
                                 c);
         Out_data->sweep[0]++;
         if (Out_data->sweep[0] >= Params.sweep_max){
            conv_flag = 1;
            break;
         }
      }
      if (conv_flag){
         break;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   Out_data->wtime_tot[0] = MPI_Wtime() - start;
   MPI_Barrier(MPI_COMM_WORLD);
}
