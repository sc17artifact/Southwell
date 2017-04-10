#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "Misc.h"
#include "MatrixUtils.h"


double RmaResNorm2(double *r, int n)
{
   double r_norm2 = Norm2(r, n);
   if (r_norm2 < DBL_EPSILON) return 0;
   else return r_norm2;
}

void DMEM_MatVecProd_CscBlocks(DMEM_MatrixInfo Mat,
                               DMEM_ProcInfo Proc, 
                               double *x,
                               double *z,
                               int n_glob)
{
   double sum, elem;
   int row, ind, neighb_rank;
   double *y = (double *)calloc(n_glob, sizeof(double));
   for (int i = 0; i < Mat.n; i++){
      for (int j = 0; j < Mat.D.j_ptr[i+1]-Mat.D.j_ptr[i]; j++){
         ind = Mat.D.j_ptr[i]+j;
         row = Mat.D.i[ind];
         elem = Mat.D.a[ind];
         y[row] += elem * x[i];
      }
   }
   for (int q = 0; q < Proc.Neighb.size; q++){
      for (int i = 0; i < Mat.n; i++){
         for (int j = 0; j < Mat.B[q].j_ptr[i+1]-Mat.B[q].j_ptr[i]; j++){
            ind = Mat.B[q].j_ptr[i]+j;
            row = Mat.B[q].i[ind];
            elem = Mat.B[q].a[ind];
            y[row] += elem * x[i];
         }
      }
   }
   MPI_Allreduce(y, z, n_glob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   free(y);
}

void DMEM_MatVecProd_MC_CscBlocks(DMEM_MulticolorInfo MC_info,
                                  double *x,
                                  double *z,
                                  int n_glob)
{
   double sum, elem;
   int row, ind, neighb_rank, disp_loc;
   double *y = (double *)calloc(n_glob, sizeof(double));
   for (int c = 0; c < MC_info.P_glob.nparts; c++){
      disp_loc = MC_info.disp_loc[c];   
      for (int i = 0; i < MC_info.D[c].n; i++){
         for (int j = 0; j < 
                 MC_info.D[c].j_ptr[i+1] - 
                    MC_info.D[c].j_ptr[i]; j++){
            ind = MC_info.D[c].j_ptr[i]+j;
            row = MC_info.D[c].i[ind];
            elem = MC_info.D[c].a[ind];
            y[row] += (elem * x[disp_loc+i]);
         }
      }
      for (int cc = 0; cc < MC_info.P_glob.nparts; cc++){
         for (int q = 0; q < MC_info.NeighbSend[c][cc].size; q++){
            for (int i = 0; i < MC_info.B[c][cc][q].n; i++){
               for (int j = 0; j < 
                       MC_info.B[c][cc][q].j_ptr[i+1] - 
                          MC_info.B[c][cc][q].j_ptr[i]; j++){
                  ind = MC_info.B[c][cc][q].j_ptr[i]+j;
                  row = MC_info.B[c][cc][q].i[ind];
                  elem = MC_info.B[c][cc][q].a[ind];
                  y[row] += (elem * x[disp_loc+i]);
               }
            }
         }
      }
   }
   MPI_Allreduce(y, z, n_glob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   free(y);
}

void DMEM_Residual_CscBlocks(DMEM_MulticolorInfo MC_info,
                             DMEM_MatrixInfo Mat, 
                             DMEM_ProcInfo Proc,
                             double *x, 
                             double *b,
                             double *r)
{
   int disp_loc;
   double *Ax;
   if (color_flag){
      Ax = (double *)calloc(MC_info.n_glob, sizeof(double));
      DMEM_MatVecProd_MC_CscBlocks(MC_info, 
                                   x, 
                                   Ax,
                                   MC_info.n_glob);
      for (int i = 0; i < MC_info.n_glob; i++){
         r[i] = b[i] - Ax[i];
      }
   }
   else {
      Ax = (double *)calloc(Mat.n_glob, sizeof(double));
      DMEM_MatVecProd_CscBlocks(Mat, Proc, x, Ax, Mat.n_glob);
      for (int i = 0; i < Mat.n_glob; i++){
         r[i] = b[i] - Ax[i];
      }
   }
   free(Ax);
}

void DMEM_Residual(DMEM_MatrixInfo Mat,
                   DMEM_MulticolorInfo MC_info,
                   DMEM_ProcInfo Proc,
                   DMEM_SolveVars *Vars)
{
   int disp_loc, disp_glob;
   int rank, world_size;
   rank = Proc.rank;
   world_size = Proc.world_size;
   DMEM_Residual_CscBlocks(MC_info,
                           Mat,
                           Proc,
                           Vars->Loc.u,
                           Vars->Glob.b,
                           Vars->Glob.r);
   if (color_flag){
      for (int c = 0; c < MC_info.P_glob.nparts; c++){
         disp_loc = MC_info.disp_loc[c];
         disp_glob = 
            MC_info.P_glob.disp[c] + MC_info.P_loc[c].disp[rank];
         for (int i = 0; i < MC_info.D[c].n; i++){
            Vars->Loc.r[disp_loc+i] = Vars->Glob.r[disp_glob+i];
         }
      }
   }
   else {
      for (int i = 0; i < Mat.n; i++){
         Vars->Loc.r[i] = Vars->Glob.r[Mat.disp+i];
      }
   }
}

void DMEM_Write_csc(char *buffer,
                    DMEM_MatrixInfo Mat,
                    DMEM_ProcInfo Proc,
                    int base)
{
   int row, col, k;
   double elem;
   int disp = Mat.disp;
   FILE *out_file;
   MPI_Barrier(MPI_COMM_WORLD);
   for (int p = 0; p < Proc.world_size; p++){
      if (p == Proc.rank){
         out_file = fopen(buffer, "a");
         for (int i = 0; i < Mat.n; i++){
            k = Mat.A.j_ptr[i];
            for (int j = 0; j < Mat.A.j_ptr[i+1]-Mat.A.j_ptr[i]; j++){
               row = Mat.A.i[k+j];
               col = i;
               elem = Mat.A.a[k+j];
               fprintf(out_file, "%d %d %e\n",
                       row+base, disp+col+base, elem);
            }
         }
         fclose(out_file);
      }
      MPI_Barrier(MPI_COMM_WORLD); 
   }
}


void DMEM_WriteBlocks_csc(char *buffer, 
                          DMEM_MatrixInfo Mat, 
                          DMEM_ProcInfo Proc, 
                          int base)
{
   int row, col, k;
   double elem;
   int disp = Mat.disp;
   FILE *out_file;
   MPI_Barrier(MPI_COMM_WORLD);
   for (int p = 0; p < Proc.world_size; p++){
      if (Proc.rank == p){
         out_file = fopen(buffer, "a");
         for (int i = 0; i < Mat.n; i++){
            k = Mat.D.j_ptr[i];
            for (int j = 0; j < Mat.D.j_ptr[i+1]-Mat.D.j_ptr[i]; j++){
               row = Mat.D.i[k+j];
               col = i;
               elem = Mat.D.a[k+j];
               fprintf(out_file, "%d %d %e\n",
                       row+base, disp+col+base, elem);
            }
         }
         for (int q = 0; q < Proc.Neighb.size; q++){
            for (int i = 0; i < Mat.n; i++){
               k = Mat.B[q].j_ptr[i];
               for (int j = 0; j < Mat.B[q].j_ptr[i+1]-Mat.B[q].j_ptr[i]; j++){
                  row = Mat.B[q].i[k+j];
                  col = i;
                  elem = Mat.B[q].a[k+j];
                  fprintf(out_file, "%d %d %e\n",
                          row+base, disp+col+base, elem);
               }
            }
         }
         fclose(out_file);
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
}

void DMEM_WriteBlocks_MC(char *buffer,
                         DMEM_MulticolorInfo MC_info,
                         DMEM_ProcInfo Proc,
                         int base)
{
   int ind, row, col, k;
   double elem;
   int disp_loc, disp_glob;
   int rank = Proc.rank;
   FILE *out_file;
   MPI_Barrier(MPI_COMM_WORLD);
   for (int p = 0; p < Proc.world_size; p++){
      if (Proc.rank == p){
         out_file = fopen(buffer, "a");
         for (int c = 0; c < MC_info.P_glob.nparts; c++){
            disp_glob = MC_info.P_glob.disp[c]; 
            disp_loc = MC_info.P_loc[c].disp[rank];
            for (int i = 0; i < MC_info.D[c].n; i++){
               for (int j = 0; j <
                       MC_info.D[c].j_ptr[i+1] -
                          MC_info.D[c].j_ptr[i]; j++){
                  col = i;
                  ind = MC_info.D[c].j_ptr[i]+j;
                  row = MC_info.D[c].i[ind];
                  elem = MC_info.D[c].a[ind];
                  fprintf(out_file, "%d %d %e\n",
                          row+base, disp_glob+disp_loc+col+base, elem);
               }
            }
            for (int cc = 0; cc < MC_info.P_glob.nparts; cc++){
               for (int q = 0; q < MC_info.NeighbSend[c][cc].size; q++){
                  for (int i = 0; i < MC_info.B[c][cc][q].n; i++){
                     for (int j = 0; j <
                             MC_info.B[c][cc][q].j_ptr[i+1] -
                                MC_info.B[c][cc][q].j_ptr[i]; j++){
                        col = i;
                        ind = MC_info.B[c][cc][q].j_ptr[i]+j;
                        row = MC_info.B[c][cc][q].i[ind];
                        elem = MC_info.B[c][cc][q].a[ind];
                        fprintf(out_file, "%d %d %e\n",
                                row+base, disp_glob+disp_loc+col+base, elem);
                     }
                  }
               }
            }
         }

         fclose(out_file);
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
}


void DMEM_DistributeVector_dbl(double *x_glob,
                               double *x,
                               DMEM_ProcInfo Proc,
                               DMEM_MatrixInfo Mat)
{
   MPI_Scatterv(x_glob,
                Mat.P.part,
                Mat.P.dispv,
                MPI_DOUBLE,
                x,
                Mat.n,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);
}

void DMEM_PrintResults(DMEM_MatrixInfo Mat,
                       DMEM_MulticolorInfo MC_info,
                       DMEM_ProcInfo Proc,
                       DMEM_SolveVars *Vars,
                       SolveData Out_data)
{
   double wtime_tot_mean; 
   double wtime_comp_mean; 
   double wtime_comm_mean;
   double wtime_conv_mean;
   double wtime_update_bound_mean;
   double wtime_comm_res_mean;
   double wtime_comm_sweep_mean;
   double wtime_res_estim_mean;
   double sweep_mean;
   double relax_scaled_tot; 
   double comm_scaled_tot; 
   double flop_scaled_tot;
   double comm_res_mean; 
   double comm_sweep_mean;
   double r_norm;
   double work_iter;
   unsigned long long sweep_tot;
   unsigned long long comm_res_tot, comm_sweep_tot;
   int relax_mask_sum;
   int n_glob;

   if (color_flag){
      n_glob = MC_info.n_glob;
   }
   else {
      n_glob = Mat.n_glob;
   }

   if (Proc.world_size == 1){
      for (int i = 0; i < Mat.n; i++){
         Vars->Glob.u[i] = Vars->Loc.u[i];
      }
   }
   else {
      if (color_flag){
         for (int c = 0; c < MC_info.P_glob.nparts; c++){
            MPI_Allgatherv(&(Vars->Loc.u[MC_info.disp_loc[c]]),
                           MC_info.D[c].n,
                           MPI_DOUBLE,
                           &(Vars->Glob.u[MC_info.P_glob.disp[c]]),
                           MC_info.P_loc[c].part,
                           MC_info.P_loc[c].disp,
                           MPI_DOUBLE,
                           MPI_COMM_WORLD);
         }
      }
      else {
         MPI_Allgatherv(Vars->Loc.u,
                        Mat.n,
                        MPI_DOUBLE,
                        Vars->Glob.u,
                        Mat.P.part,
                        Mat.P.dispv,
                        MPI_DOUBLE,
                        MPI_COMM_WORLD);
      }
   }
   DMEM_Residual(Mat, MC_info, Proc, Vars);
   MPI_Reduce(&(Out_data.relax_mask[0]), &relax_mask_sum, 1,
              MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.sweep[0]), &sweep_tot, 1, 
              MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.relax_scaled[0]), &relax_scaled_tot, 1, 
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.comm_scaled[0]), &comm_scaled_tot, 1, 
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.flop_scaled[0]), &flop_scaled_tot, 1,
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.wtime_tot[0]), &wtime_tot_mean, 1, 
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.wtime_comp[0]), &wtime_comp_mean, 1,
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.wtime_comm[0]), &wtime_comm_mean, 1,
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.wtime_conv[0]), &wtime_conv_mean, 1,
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&(Out_data.wtime_update_bound[0]), &wtime_update_bound_mean, 1,
              MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   if (solver_flag == DMEM_SOLVER_DS ||
       solver_flag == DMEM_SOLVER_PS){
      MPI_Reduce(&(Out_data.wtime_comm_res[0]), &wtime_comm_res_mean, 1,
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Out_data.wtime_comm_sweep[0]), &wtime_comm_sweep_mean, 1,
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Out_data.comm_res[0]), &comm_res_tot, 1, 
                 MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&(Out_data.comm_sweep[0]), &comm_sweep_tot, 1, 
                 MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   }
   if (solver_flag == DMEM_SOLVER_DS){
      MPI_Reduce(&(Out_data.wtime_res_estim[0]), &wtime_res_estim_mean, 1,
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   }
   if (Proc.rank == 0){
     // for (int i = 0; i < Mat.n_glob; i++)
     //    printf("%e\n", Vars->Glob.u[i]);
      work_iter = (double)relax_mask_sum/(double)Proc.world_size;
      sweep_mean = (double)sweep_tot/(double)Proc.world_size;
      comm_sweep_mean = (double)comm_sweep_tot/(double)Proc.world_size;
      comm_res_mean = (double)comm_res_tot/(double)Proc.world_size;
      wtime_tot_mean /= (double)Proc.world_size;
      wtime_comp_mean /= (double)Proc.world_size;
      wtime_comm_mean /= (double)Proc.world_size;
      wtime_conv_mean /= (double)Proc.world_size;
      wtime_comm_res_mean /= (double)Proc.world_size;
      wtime_comm_sweep_mean /= (double)Proc.world_size;
      wtime_update_bound_mean /= (double)Proc.world_size;
      wtime_res_estim_mean /= (double)Proc.world_size;
      r_norm = Norm2(Vars->Glob.r, n_glob);
      
      if (format_out_flag){
         printf("%e %e %e %e %e %e %e %e %e %e ", 
                r_norm,
                wtime_tot_mean,
                wtime_comp_mean,
                wtime_comm_mean,
                wtime_conv_mean,
                sweep_mean,
                comm_scaled_tot,
                relax_scaled_tot,
                work_iter,
                wtime_update_bound_mean);
         if (solver_flag == DMEM_SOLVER_DS ||
             solver_flag == DMEM_SOLVER_PS){
            printf("%e %e %e %e ",
                   comm_sweep_mean,
                   comm_res_mean,
                   wtime_comm_sweep_mean,
                   wtime_comm_res_mean);
         }
         if (solver_flag == DMEM_SOLVER_DS){
            printf("%e",
                   wtime_res_estim_mean);
         }
         printf("\n");
      }
      else {
         printf("\n\tresidual 2-norm = %e\n", r_norm);
         printf("\n\tTimings:\n\t\ttotal wall-time = %e"
                "\n\t\tcomputation time/p = %e\n\t\tcommunication time/p = %e\n\t\tmean convergence test time = %e\n\t\tupdating boundary points time/p = %e",
                wtime_tot_mean,
                wtime_comp_mean,
                wtime_comm_mean,
                wtime_conv_mean,
                wtime_update_bound_mean);
         if (solver_flag == DMEM_SOLVER_DS ||
             solver_flag == DMEM_SOLVER_PS){
            printf("\n\t\tsweep comm time/p = %e\n\t\tres comm time/p = %e",
                   wtime_comm_sweep_mean,
                   wtime_comm_res_mean);
         }
         if (solver_flag == DMEM_SOLVER_DS){
            printf("\n\t\tres estim time/p = %e",
                   wtime_res_estim_mean);
         }
         printf("\n\n\tOperations:\n\t\tsweeps = %e\n\t\ttotal comm/p = %e\n\t\trelaxations/n = %e\n\t\tfraction of active processes = %e",
                sweep_mean,
                comm_scaled_tot,
                relax_scaled_tot,
                work_iter);
         if (solver_flag == DMEM_SOLVER_DS ||
             solver_flag == DMEM_SOLVER_PS){
            printf("\n\t\tsweep comm/p = %e\n\t\tres comm/p = %e",
                   comm_sweep_mean,
                   comm_res_mean);
         }
         printf("\n");
      }
   }
}

double DMEM_BoundSumSqu(DMEM_ProcInfo *Proc,
                        DMEM_SolveVars Vars)
{
   int point;
   double sum_squ = 0;
   for (int q = 0; q < Proc->Neighb.size; q++){
      for (int i = 0; i < Proc->Bound.num_recv_points[q]; i++){
         point = Proc->Bound.recv_points[q][i];
         sum_squ += pow(fabs(Vars.Loc.r[point]), 2);
      }
   } 
   return sum_squ;
}

void DMEM_RandDouble(DMEM_ProcInfo Proc,
                     double *v,
                     int n,
                     double low,
                     double high)
{
   VSLStreamStatePtr stream;
   vslNewStream(&stream, VSL_BRNG_SFMT19937, Proc.rank);
   vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, v, -.5, .5);
}

void DMEM_Free_SolveVars(DMEM_SolveVars *Vars)
{
   free(Vars->Loc.r);
   free(Vars->Loc.u);
   free(Vars->Loc.b);
   free(Vars->Loc.du);
   free(Vars->Glob.r);
   free(Vars->Glob.u);
   free(Vars->Glob.b);
}

void DMEM_Free_MatrixInfo(DMEM_MatrixInfo *Mat,
                          int neighb_size)
{
   FreeCSC(&(Mat->D));
   for (int q = 0; q < neighb_size; q++){
      free(Mat->B[q].j_ptr);
      free(Mat->B[q].i);
      free(Mat->B[q].a);
   }
   free(Mat->B);
   FreeOrdering(&(Mat->P));
}

void DMEM_Free_ConvInfo(DMEM_ConvergeInfo *Conv)
{
   MPI_Win_free(&(Conv->win));
}

void DMEM_Free_ResNormInfo(DMEM_ResNormInfo *Res)
{
   free(Res->norm);
   if (solver_flag == DMEM_SOLVER_DS){
      free(Res->neighb_norm_estim);
      free(Res->norm_estim_squ);
      free(Res->my_norm_estim_squ);
      free(Res->my_norm_estim_squ_prev);
   }
   else if ((solver_flag == DMEM_SOLVER_PS) && (implement_flag == DMEM_POS)){
      free(Res->norm_prev);
   }
   if ((solver_flag == DMEM_SOLVER_PS) && (implement_flag == DMEM_SOS)){
      free(Res->Rma.send_disp);
      MPI_Win_free(&(Res->Rma.win));
   }
}

void DMEM_Free_BoundInfo(DMEM_BoundInfo *Bound,
                         int neighb_size)
{
   free(Bound->num_send_points);
   free(Bound->row_min);
   free(Bound->row_max);
   free(Bound->row_range_minmax);
   free(Bound->row_range_disp);
   free(Bound->Rma.send_disp);
   free(Bound->Rma.recv_disp);
   free(Bound->Rma.send_count);
   for (int q = 0; q < neighb_size; q++){
      free(Bound->points[q]);
      free(Bound->send_points[q]);
      free(Bound->recv_points[q]);
      free(Bound->Rma.send[q]);
      if (solver_flag == DMEM_SOLVER_DS){
         free(Bound->res[q]);
         free(Bound->Rma.send_prev[q]);
      }
   }
   free(Bound->send_points);
   free(Bound->recv_points);
   free(Bound->points);
   free(Bound->Rma.send);
   if (solver_flag == DMEM_SOLVER_DS){
      free(Bound->res);
      free(Bound->Rma.send_prev);
   }
   MPI_Win_free(&(Bound->Rma.win));
}

void DMEM_Free_NeighbInfo(DMEM_NeighbInfo *Neighb)
{
   free(Neighb->ranks);
}

void DMEM_Free_ProcInfo(DMEM_ProcInfo *Proc)
{
   DMEM_Free_ResNormInfo(&(Proc->Res));
   if (implement_flag == DMEM_POS){
      DMEM_Free_ConvInfo(&Proc->Conv);
   }
   DMEM_Free_BoundInfo(&Proc->Bound, Proc->Neighb.size);
   DMEM_Free_NeighbInfo(&Proc->Neighb);
}
