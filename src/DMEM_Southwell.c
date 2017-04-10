#include "Southwell.h"
#include "DMEM_Southwell.h"
#include "SEQ_Solve.h"
#include "DMEM_Solve.h"
#include "DMEM_Setup.h"
#include "DMEM_Misc.h"
#include "Misc.h"
#include "FileUtils.h"
#include "DMEM_MCGSutils.h"

int mat_file_flag = 0;
int format_out_flag = 0;
int color_flag = 0;
int reorder_flag = 1;
int solver_flag = DMEM_SOLVER_J;
int loc_solver_flag = DMEM_LOC_SOLVER_SEQGS;
int implement_flag = DMEM_POS;
int ds_res_estim_flag = 1;
int no_solver = 1;
int num_samps = 1;

int main(int argc, char *argv[])
{
   double low = -1, high = 1, rand_range = high-low;
   double rand_num;
   double r_norm;
   double x_low = -1, x_high = 1, b_low = -1, b_high = 1;
   double x_sin_scale, b_sin_scale;
   double b_norm, Ax_norm;
   double tol_high, tol_low, tol_incr;
   double start;
   int sweep_low, sweep_high, sweep_incr;
   int r_norm_order; 
   int disp_loc, disp_glob;
   int n, n_glob;
 
   int delay_relax = 0;
   int r_unit_norm_flag = 1;
   int x_file_flag = 0, b_file_flag = 0;
   int x_zeros_flag = 0, b_zeros_flag = 0;
   int x_ones_flag = 0, b_ones_flag = 0;
   int x_sin_flag = 0, b_sin_flag = 0;
   int solver_seq_s = 0; 
   int solver_seq_gs = 0;
   int solver_sos_sj = 0, solver_pos_aj = 0;
   int solver_sos_sps = 0, solver_pos_aps = 0;
   int solver_sos_sds = 0, solver_pos_ads = 0;
   int solver_sos_mcgs = 0;
   int comm_out_flag = 0, time_out_flag = 0;
   
   MetisGraph G;
   DMEM_MatrixInfo *Mat = 
      (DMEM_MatrixInfo *)malloc(sizeof(DMEM_MatrixInfo));
   DMEM_MulticolorInfo *MC_info = 
      (DMEM_MulticolorInfo *)malloc(sizeof(DMEM_MulticolorInfo));
   DMEM_ProcInfo Proc;
   SolveParams Params;
   DMEM_SolveVars Vars;
   SolveData Out_data;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &(Proc.world_size));
   MPI_Comm_rank(MPI_COMM_WORLD, &(Proc.rank));
   omp_set_dynamic(0);
   omp_set_num_threads(1);

   Params.tol = tol_low = 1e-15;
   Params.sweep_max = 0;
   Params.ds_delay = 0;
   int arg_iter = 0;
   int m1 = 100, m2 = 100;
   
   char mat_str[100] = "";
   char x_file_str[100] = "";
   char b_file_str[100] = "";
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-solver") == 0){
         arg_iter++;
         if (strcmp(argv[arg_iter], "sos_sj") == 0){
            no_solver = 0;
            solver_sos_sj = 1;
            solver_flag = DMEM_SOLVER_J;
            implement_flag = DMEM_SOS;
         }
         else if (strcmp(argv[arg_iter], "pos_aj") == 0){
            no_solver = 0;
            solver_pos_aj = 1;
            solver_flag = DMEM_SOLVER_J;
            implement_flag = DMEM_POS;
         }
         else if (strcmp(argv[arg_iter], "sos_sps") == 0){
            no_solver = 0;
            solver_sos_sps = 1;
            solver_flag = DMEM_SOLVER_PS;
            implement_flag = DMEM_SOS;
         }
         else if (strcmp(argv[arg_iter], "pos_aps") == 0){
            no_solver = 0;
            solver_pos_aps = 1;
            solver_flag = DMEM_SOLVER_PS;
            implement_flag = DMEM_POS;
         }
         else if (strcmp(argv[arg_iter], "sos_sds") == 0){
            no_solver = 0;
            solver_sos_sds = 1;
            solver_flag = DMEM_SOLVER_DS;
            implement_flag = DMEM_SOS;
         }
         else if (strcmp(argv[arg_iter], "pos_ads") == 0){
            no_solver = 0;
            solver_pos_ads = 1;
            solver_flag = DMEM_SOLVER_DS;
            implement_flag = DMEM_POS;
         }
         else if (strcmp(argv[arg_iter], "sos_mcgs") == 0){
            no_solver = 0;
            solver_sos_mcgs = 1;
            color_flag = 1;
            solver_flag = DMEM_SOLVER_MCGS;
            implement_flag = DMEM_SOS;
         }
      }
      if (strcmp(argv[arg_iter], "-loc_solver") == 0){
         arg_iter++;
         if (strcmp(argv[arg_iter], "direct") == 0){
            loc_solver_flag = DMEM_LOC_SOLVER_DIRECT;
         }
         else if (strcmp(argv[arg_iter], "seq_gs") == 0){
            loc_solver_flag = DMEM_LOC_SOLVER_SEQGS;
         }
      }
      else if (strcmp(argv[arg_iter], "-tol") == 0){
         arg_iter++;
         Params.tol = tol_low = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_max") == 0){
         arg_iter++;
         Params.sweep_max = sweep_low = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-num_samps") == 0){
         arg_iter++;
         num_samps = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-ds_delay") == 0){
         arg_iter++;
         Params.ds_delay = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-ds_no_res_estim") == 0){
         ds_res_estim_flag = 0;
      }
      else if (strcmp(argv[arg_iter], "-mat_file") == 0){
         arg_iter++;
         strcpy(mat_str, argv[arg_iter]);
         mat_file_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-r_unit_norm") == 0){
         r_unit_norm_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-x_file") == 0){
         arg_iter++;
         strcpy(x_file_str, argv[arg_iter]);
         x_file_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-x_zeros") == 0){
          x_zeros_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-x_ones") == 0){
          x_ones_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-x_sin") == 0){
          x_sin_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-x_rand_range") == 0){
          arg_iter++;
          x_low = atof(argv[arg_iter++]);
          x_high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-b_file") == 0){
         arg_iter++;
         strcpy(b_file_str, argv[arg_iter]);
         b_file_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-b_zeros") == 0){
          b_zeros_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-b_ones") == 0){
          b_ones_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-b_sin") == 0){
          b_sin_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-b_rand_range") == 0){
          arg_iter++;
          b_low = atof(argv[arg_iter++]);
          b_high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-Laplace2D_FD5pt_size") == 0){
         arg_iter++;
         m1 = atoi(argv[arg_iter++]);
         m2 = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-format_out") == 0){
         format_out_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-no_reorder") == 0){
         reorder_flag = 0;
      }
      arg_iter++;
   }

   tol_incr = 100;
   tol_high = tol_low/2;
   arg_iter = 0;
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-tol_low") == 0){
         arg_iter++;
         tol_low = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-tol_high") == 0){
         arg_iter++;
         tol_high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-tol_incr") == 0){
         arg_iter++;
         tol_incr = atof(argv[arg_iter]);
      }
      arg_iter++; 
   } 

   arg_iter = 0;
   sweep_high = sweep_low; 
   sweep_incr = (int)1e9;
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-sweep_low") == 0){
         arg_iter++;
         sweep_low = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_high") == 0){
         arg_iter++;
         sweep_high = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_incr") == 0){
         arg_iter++;
         sweep_incr = atoi(argv[arg_iter]);
      }
      arg_iter++;
   }

   FILE *mat_file = fopen(mat_str, "rb");
   start = MPI_Wtime();

   if (color_flag){
      DMEM_MulticolorSetup(mat_file, MC_info, &Proc, m1, m2);
   }
   else {
      DMEM_MetisSetup(mat_file, Mat, &Proc, m1, m2);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if ((Proc.rank == 0) && !format_out_flag){
      printf("~~~~~~~~~~~~~~~~~~~~~~~~\n"
             "~ TOTAL SETUP TIME = %e \n" 
             "~~~~~~~~~~~~~~~~~~~~~~~~\n", 
             MPI_Wtime() - start);
   }
   SetupOutData(&Out_data, 0);
   x_sin_scale = b_sin_scale = Mat->n_glob;
   if (color_flag){
      n = MC_info->n;
      n_glob = MC_info->n_glob;
   }
   else {
      n = Mat->n;
      n_glob = Mat->n_glob;
   }
   Vars.Loc.r = (double *)calloc(n, sizeof(double));
   Vars.Loc.b = (double *)calloc(n, sizeof(double));
   Vars.Loc.u = (double *)calloc(n, sizeof(double));
   Vars.Loc.du = (double *)calloc(n, sizeof(double));
   Vars.Glob.u = (double *)calloc(n_glob, sizeof(double));
   Vars.Glob.r = (double *)calloc(n_glob, sizeof(double));
   Vars.Glob.b = (double *)calloc(n_glob, sizeof(double));
   double *u_init = (double *)calloc(n_glob, sizeof(double));
   double *r_init = (double *)calloc(n_glob, sizeof(double));
   
   double *ones = (double *)calloc(n, sizeof(double));
   for (int i = 0; i < n; i++) ones[i] = 1.0;
 
   if (x_file_flag){
      if (Proc.rank == 0){
         ReadVector(x_file_str, &Vars.Glob.u, n_glob);
      }
      DMEM_DistributeVector_dbl(Vars.Glob.u, Vars.Loc.u, Proc, *Mat);
   }
   else if (x_zeros_flag){
      for (int i = 0; i < n; i++) Vars.Loc.u[i] = 0;
   }
   else if (x_ones_flag){
      for (int i = 0; i < n; i++) Vars.Loc.u[i] = 1.0;
   }
   else if (x_sin_flag){
      for (int i = 0; i < n; i++) 
         Vars.Loc.u[i] = x_sin_scale * sin((double)i);
   }
   else {
      DMEM_RandDouble(Proc, Vars.Loc.u, n, x_low, x_high);
   }
 
   if (b_file_flag){
      if (Proc.rank == 0){
         ReadVector(b_file_str, &Vars.Glob.b, n_glob);
      }
      DMEM_DistributeVector_dbl(Vars.Glob.b, Vars.Loc.b, Proc, *Mat);
   }
   else if (b_zeros_flag){
      for (int i = 0; i < n; i++) Vars.Loc.b[i] = 0;
   }
   else if (b_ones_flag){
      for (int i = 0; i < n; i++) Vars.Loc.b[i] = 1.0;
   }
   else if (b_sin_flag){
      for (int i = 0; i < n; i++)
         Vars.Loc.b[i] = b_sin_scale * sin((double)i);
   }
   else {
      DMEM_RandDouble(Proc, Vars.Loc.b, n, x_low, x_high);
   }
 
   MPI_Barrier(MPI_COMM_WORLD);
   if (color_flag){
      for (int c = 0; c < MC_info->P_glob.nparts; c++){
         MPI_Allgatherv(&(Vars.Loc.u[MC_info->disp_loc[c]]),
                        MC_info->D[c].n,
                        MPI_DOUBLE,
                        &(Vars.Glob.u[MC_info->P_glob.disp[c]]),
                        MC_info->P_loc[c].part,
                        MC_info->P_loc[c].disp,
                        MPI_DOUBLE,
                        MPI_COMM_WORLD);
         MPI_Allgatherv(&(Vars.Loc.b[MC_info->disp_loc[c]]),
                        MC_info->D[c].n,
                        MPI_DOUBLE,
                        &(Vars.Glob.b[MC_info->P_glob.disp[c]]),
                        MC_info->P_loc[c].part,
                        MC_info->P_loc[c].disp,
                        MPI_DOUBLE,
                        MPI_COMM_WORLD);
      }
   }
   else {
      MPI_Allgatherv(Vars.Loc.u,
                     Mat->n,
                     MPI_DOUBLE,
                     Vars.Glob.u,
                     Mat->P.part,
                     Mat->P.dispv,
                     MPI_DOUBLE,
                     MPI_COMM_WORLD);
      MPI_Allgatherv(Vars.Loc.b,
                     Mat->n,
                     MPI_DOUBLE,
                     Vars.Glob.b,
                     Mat->P.part,
                     Mat->P.dispv,
                     MPI_DOUBLE,
                     MPI_COMM_WORLD);
   }
   DMEM_Residual(*Mat, *MC_info, Proc, &Vars);
 
   start = MPI_Wtime();
   if (r_unit_norm_flag){
      if (b_zeros_flag){
         if (color_flag){
            double *Ax = (double *)calloc(n_glob, sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            DMEM_MatVecProd_MC_CscBlocks(*MC_info,
                                         Vars.Loc.u,
                                         Ax,
                                         n_glob);
            MPI_Barrier(MPI_COMM_WORLD);
            Ax_norm = Norm2(Ax, n_glob);
            free(Ax);
            for (int i = 0; i < n; i++)
               Vars.Loc.u[i] /= Ax_norm;
            MPI_Barrier(MPI_COMM_WORLD);
            for (int c = 0; c < MC_info->P_glob.nparts; c++){
               disp_loc = MC_info->disp_loc[c];
               disp_glob = MC_info->P_glob.disp[c];
               MPI_Allgatherv(&(Vars.Loc.u[MC_info->disp_loc[c]]),
                              MC_info->D[c].n,
                              MPI_DOUBLE,
                              &(Vars.Glob.u[MC_info->P_glob.disp[c]]),
                              MC_info->P_loc[c].part,
                              MC_info->P_loc[c].disp,
                              MPI_DOUBLE,
                              MPI_COMM_WORLD);
            }
         }
         else {
            double *Ax = (double *)calloc(n_glob, sizeof(double));
            int count = 1, mat_power = 2;
            DMEM_MatVecProd_CscBlocks(*Mat, Proc, Vars.Loc.u, Ax, n_glob);
            while (count < mat_power){
               for (int i = 0; i < n; i++)
                  Vars.Loc.u[i] = Ax[Mat->disp+i];
               DMEM_MatVecProd_CscBlocks(*Mat, Proc, Vars.Loc.u, Ax, n_glob);
               count++;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            Ax_norm = Norm2(Ax, n_glob);
            free(Ax);
            for (int i = 0; i < n; i++)
               Vars.Loc.u[i] /= Ax_norm;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allgatherv(Vars.Loc.u,
                           Mat->n,
                           MPI_DOUBLE,
                           Vars.Glob.u,
                           Mat->P.part,
                           Mat->P.dispv,
                           MPI_DOUBLE,
                           MPI_COMM_WORLD);
         }
         MPI_Barrier(MPI_COMM_WORLD);
         DMEM_Residual(*Mat, *MC_info, Proc, &Vars);
         MPI_Barrier(MPI_COMM_WORLD);
      }
      else if (x_zeros_flag){
         b_norm = Norm2(Vars.Glob.b, n_glob);
         for (int i = 0; i < n; i++)
            Vars.Loc.b[i] /= b_norm;
         if (color_flag){
            for (int c = 0; c < MC_info->P_glob.nparts; c++){
               MPI_Allgatherv(&(Vars.Loc.b[MC_info->disp_loc[c]]),
                              MC_info->D[c].n,
                              MPI_DOUBLE,
                              &(Vars.Glob.b[MC_info->P_glob.disp[c]]),
                              MC_info->P_loc[c].part,
                              MC_info->P_loc[c].disp,
                              MPI_DOUBLE,
                              MPI_COMM_WORLD);
            }
            
         }
         else {
            MPI_Allgatherv(Vars.Loc.b,
                           Mat->n,
                           MPI_DOUBLE,
                           Vars.Glob.b,
                           Mat->P.part,
                           Mat->P.dispv,
                           MPI_DOUBLE,
                           MPI_COMM_WORLD);
         }
         DMEM_Residual(*Mat, *MC_info, Proc, &Vars);
      }
      if ((Proc.rank == 0) && !format_out_flag){
         printf("\nInitial guess scaled, "
                "initial residual norm = %e, time = %e\n",
                Norm2(Vars.Glob.r, n_glob), MPI_Wtime() - start);
      }
   }
 // if (color_flag){
 //    DMEM_MatVecProd_MC_CscBlocks(*MC_info,
 //                                 ones,
 //                                 Vars.Glob.b,
 //                                 n_glob);
 // }
 // else {
 //    DMEM_MatVecProd_CscBlocks(*Mat,
 //                              Proc,
 //                              ones,
 //                              Vars.Glob.b,
 //                              n_glob);
 // }
  DMEM_Residual(*Mat, *MC_info, Proc, &Vars);
  for (int i = 0; i < n_glob; i++){
     u_init[i] = Vars.Glob.u[i];
     r_init[i] = Vars.Glob.r[i];
  }
  Proc.res_norm_init = Norm2(Vars.Glob.r, n_glob);

   for (int sweep = sweep_low; sweep <= sweep_high; sweep += sweep_incr){
      Params.sweep_max = sweep;
      for (int samp = 0; samp < num_samps; samp++){
         for (int i = 0; i < Mat->n; i++){
            Vars.Loc.r[i] = r_init[Mat->disp+i];
            Vars.Loc.u[i] = u_init[Mat->disp+i];
         }
         if (color_flag){
            for (int c = 0; c < MC_info->P_glob.nparts; c++){
               DMEM_MC_ZeroData(MC_info);
            }
         }
         else {
            DMEM_ZeroProc(&Proc);
         } 
         ZeroOutData(&Out_data, 0);
         if (loc_solver_flag == DMEM_LOC_SOLVER_DIRECT){
            Out_data.flop_scaled[0] =
               (double)(Mat->Pard.iparm[18]*1e6)/(double)Mat->n_glob;
         }
         MPI_Barrier(MPI_COMM_WORLD);
         if (solver_pos_aj){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nASYNCHRONOUS JACOBI USING PASSIVE ONE-SIDED MPI\n");
            }
            DMEM_POS_AsyncBlockJacobi(*Mat, Params, &Proc, &Vars, &Out_data); 
         }
         if (solver_sos_sj){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nSYNCHRONOUS JACOBI USING SCALABLE ONE-SIDED MPI\n");
            }
            DMEM_SOS_SyncBlockJacobi(*Mat, Params, &Proc, &Vars, &Out_data);
         }
         if (solver_pos_ads){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nASYNCHRONOUS DISTRIBUTED SOUTHWELL USING PASSIVE ONE-SIDED MPI\n");
            }
            DMEM_POS_AsyncDistrSouthwell(*Mat, Params, &Proc, &Vars, &Out_data);
         }
         if (solver_sos_sds){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nSYNCHRONOUS DISTRIBUTED SOUTHWELL USING SCALABLE ONE-SIDED MPI\n");
            }
            DMEM_SOS_SyncDistrSouthwell(*Mat, Params, &Proc, &Vars, &Out_data);
         }
         if (solver_sos_sps){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nSYNCHRONOUS PARALLEL SOUTHWELL USING SCALABLE ONE-SIDED MPI\n");
            }
            DMEM_SOS_SyncParSouthwell(*Mat, Params, &Proc, &Vars, &Out_data);
         }
         if (solver_pos_aps){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nASYNCHRONOUS PARALLEL SOUTHWELL USING PASSIVE ONE-SIDED MPI\n");
            }
            DMEM_POS_AsyncParSouthwell(*Mat, Params, &Proc, &Vars, &Out_data);
         }
         if (solver_sos_mcgs){
            if ((Proc.rank == 0) && !format_out_flag){
               printf("\nSYNCHRONOUS MCGS USING SCALABLE ONE-SIDED MPI\n");
            }
            DMEM_SOS_MulticolorGS(MC_info, Proc, Params, &Vars, &Out_data);
         }
         if (!no_solver){
            DMEM_PrintResults(*Mat, *MC_info, Proc, &Vars, Out_data);
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   if (color_flag){
   }
   else {
      DMEM_Free_ProcInfo(&Proc);
      DMEM_Free_MatrixInfo(Mat, Proc.Neighb.size);
   }
   
   MPI_Finalize();
   return 0;
}
