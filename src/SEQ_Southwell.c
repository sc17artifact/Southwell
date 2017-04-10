#include "Southwell.h"
#include "SEQ_Solve.h"
#include "SEQ_Setup.h"
#include "Misc.h"
#include "FileUtils.h"

int mat_file_flag = 0;
int format_out_flag = 0;
int color_flag = 0;
int threads = 1;
int block_flag = 0;
int loc_direct_flag = 0;
int ps_num_relax = 1;

int main(int argc, char *argv[])
{
   double low = -1, high = 1, rand_range = high-low;
   double rand_num;
   double r_norm;
   double x_low = -1, x_high = 1, b_low = -1, b_high = 1;
   double x_sin_scale, b_sin_scale;
   double b_norm, Ax_norm;
   double tol_high, tol_low, tol_incr_mult;
   double start;

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
   int solver_seq_mcgs = 0;
   int solver_seq_ps = 0;
   int solver_seq_j = 0;
   int solver_seq_bj = 0;
   int solver_seq_bps = 0;
   int samp_low, samp_high;
   int sweep_low, sweep_high, sweep_incr_add; 
   
  
   PardisoInfo *Pard;
   CSC A, *B, *D;
   OrderInfo P;
   P.nparts = 2;
   SolveParams Params;
   SolveVars Vars;
   SolveData Out_data;

   samp_low = 0; samp_high = 1;
   Params.tol = tol_low = 1e-3;
   Params.sweep_max = sweep_low = (int)1e9;
   Params.ds_delay = 0;
   int arg_iter = 0;
   int nx = 100, ny = 100;
   
   char buffer[256];
   char mat_str[100] = "";
   char x_file_str[100] = "";
   char b_file_str[100] = "";
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-solver") == 0){
         arg_iter++;
         if (strcmp(argv[arg_iter], "s") == 0){
            solver_seq_s = 1;
         }
         else if (strcmp(argv[arg_iter], "gs") == 0){
            solver_seq_gs = 1;
         }
         else if (strcmp(argv[arg_iter], "mcgs") == 0){
            solver_seq_mcgs = 1;
            color_flag = 1;
         }
         else if (strcmp(argv[arg_iter], "ps") == 0){
            solver_seq_ps = 1;
         }
         else if (strcmp(argv[arg_iter], "j") == 0){
            solver_seq_j = 1;
         }
         else if (strcmp(argv[arg_iter], "bj") == 0){
            block_flag = 1;
            solver_seq_bj = 1;
         }
         else if (strcmp(argv[arg_iter], "bps") == 0){
            block_flag = 1;
            solver_seq_bps = 1;
         }
      }
      else if (strcmp(argv[arg_iter], "-loc_solver") == 0){
         arg_iter++;
         if (strcmp(argv[arg_iter], "direct") == 0){
            loc_direct_flag = 1;
         }
      }
      else if (strcmp(argv[arg_iter], "-ps_num_relax") == 0){
         arg_iter++;
         ps_num_relax = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-tol") == 0){
         arg_iter++;
         Params.tol = tol_low = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_max") == 0){
         arg_iter++;
         Params.sweep_max = sweep_low = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-num_blocks") == 0){
         arg_iter++;
         P.nparts = atoi(argv[arg_iter]);
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
         nx = atoi(argv[arg_iter++]);
         ny = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-format_out") == 0){
         format_out_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-samp_low") == 0){
         arg_iter++;
         samp_low = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-samp_high") == 0){
         arg_iter++;
         samp_high = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-threads") == 0){
         arg_iter++;
         threads = atoi(argv[arg_iter]);
      }
      arg_iter++;
   }

   tol_incr_mult = 100;
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
      else if (strcmp(argv[arg_iter], "-tol_incr_mult") == 0){
         arg_iter++;
         tol_incr_mult = atof(argv[arg_iter]);
      }
      arg_iter++; 
   } 

   sweep_incr_add = (int)1e9;
   sweep_high = sweep_low;
   arg_iter = 0;
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-sweep_low") == 0){
         arg_iter++;
         sweep_low = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_high") == 0){
         arg_iter++;
         sweep_high = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_incr_add") == 0){
         arg_iter++;
         sweep_incr_add = atoi(argv[arg_iter]);
      }
      arg_iter++;
   }

   FILE *mat_file = fopen(mat_str, "rb");
   start = omp_get_wtime();
   
   SEQ_Setup(mat_file, &A, &B, &D, &P, &Pard, nx, ny);
  // strcpy(buffer, "metis_matrix_matlab.txt");
  // WriteBlocks_csc(buffer, D, B, P, 1);

   if (!format_out_flag){
      printf("~~~~~~~~~~~~~~~~~~~~~~~~\n"
             "~ TOTAL SETUP TIME = %e \n" 
             "~~~~~~~~~~~~~~~~~~~~~~~~\n", 
             omp_get_wtime() - start);
   }
   n = A.n;
   SetupOutData(&Out_data, 1);

   for (int sweep = sweep_low; sweep <= sweep_high; sweep += sweep_incr_add){
      Params.sweep_max = sweep;
  // for (double tol = tol_low; tol > tol_high; tol /= tol_incr_mult){
  //    Params.tol = tol;
      ZeroOutData(&Out_data, 0);         
   
      x_sin_scale = b_sin_scale = n;
      Vars.r = (double *)calloc(n, sizeof(double));
      Vars.b = (double *)calloc(n, sizeof(double));
      Vars.x = (double *)calloc(n, sizeof(double));
      Vars.u = (double *)calloc(n, sizeof(double));
      Vars.du = (double *)calloc(n, sizeof(double));
   
      if (x_file_flag){
         ReadVector(x_file_str, &Vars.u, n);
      }
      else if (x_zeros_flag){
         for (int i = 0; i < n; i++) Vars.u[i] = 0;
      }
      else if (x_ones_flag){
         for (int i = 0; i < n; i++) Vars.u[i] = 1.0;
      }
      else if (x_sin_flag){
         for (int i = 0; i < n; i++) 
            Vars.u[i] = x_sin_scale * sin((double)i);
      }
      else {
         RandDouble(Vars.u, n, x_low, x_high);
      }
   
      if (b_file_flag){
         ReadVector(b_file_str, &Vars.b, n);
      }
      else if (b_zeros_flag){
         for (int i = 0; i < n; i++) Vars.b[i] = 0;
      }
      else if (b_ones_flag){
         for (int i = 0; i < n; i++) Vars.b[i] = 1.0;
      }
      else if (b_sin_flag){
         for (int i = 0; i < n; i++)
            Vars.b[i] = b_sin_scale * sin((double)i);
      }
      else {
         RandDouble(Vars.b, n, x_low, x_high);
      }
      Residual(A, Vars.u, Vars.b, Vars.r);
      RandDouble(Vars.x, n, x_low, x_high);
      double *ones = (double *)calloc(n, sizeof(double));
      for (int i = 0; i < n; i++) ones[i] = 1;
      MatVecProd(A, Vars.x, Vars.b);
    
      start = omp_get_wtime();
      if (r_unit_norm_flag){
         if (b_zeros_flag){
            int count = 1, mat_power = 1;
            double *Ax = (double *)calloc(n, sizeof(double));
            MatVecProd(A, Vars.u, Ax);
            while (count < mat_power){
               for (int i = 0; i < n; i++)
                  Vars.u[i] = Ax[i];
               MatVecProd(A, Vars.u, Ax);
               count++;
            }
            MatVecProd(A, Vars.u, Ax);
            Ax_norm = Norm2(Ax, n);
            free(Ax);
            for (int i = 0; i < n; i++)
               Vars.u[i] /= Ax_norm;
         }
         else if (x_zeros_flag){
            b_norm = Norm2(Vars.b, n);
            for (int i = 0; i < n; i++)
               Vars.b[i] /= b_norm;
         }
         Residual(A, Vars.u, Vars.b, Vars.r);
         if (!format_out_flag){
            printf("\nInitial guess scaled, "
                   "initial residual norm = %e, time = %e\n",
                   Norm2(Vars.r, n), omp_get_wtime() - start);
         }
      }
      
      for (int samp = samp_low; samp < samp_high; samp++){
         if (solver_seq_s){
            Out_data.relax_hist = (int *)calloc(n, sizeof(int));
            if (!format_out_flag){
               printf("\nSEQUENTIAL SOUTHWELL\n");
            }
            SEQ_Southwell(A, &Vars, Params, &Out_data); 
            free(Out_data.relax_hist);
         }
         if (solver_seq_gs){
            if (!format_out_flag){
               printf("\nGAUSS-SEIDEL\n");
            }
            SEQ_GaussSeidelRes(A, &Vars, Params, &Out_data);
         }
         if (solver_seq_mcgs){
            if (!format_out_flag){
               printf("\nMULTICOLOR GAUSS-SEIDEL\n");
            }
            SMEM_MCGS(A, P, &Vars, Params, &Out_data);
         }
         if (solver_seq_ps){
            int *degree = (int *)calloc(n, sizeof(int));
            DegreeCSC(A, degree);
            Out_data.relax_hist = (int *)calloc(n, sizeof(int));
            Out_data.relax_mask = (int *)calloc(n, sizeof(int));
            if (!format_out_flag){
               printf("\nPARALLEL SOUTHWELL\n");
            }
            SMEM_ParSouthwellRes(A, &Vars, Params, &Out_data);
           // for (int i = 0; i < n; i++){
           //    printf("%d, %d\n", degree[i], Out_data.relax_hist[i]); 
           // }
            free(Out_data.relax_hist);
            free(Out_data.relax_mask);
            free(degree);
         }
         if (solver_seq_j){
            if (!format_out_flag){
               printf("\nJACOBI\n");
            }
            SMEM_JacobiRes(A, &Vars, Params, &Out_data);
         }
         if (solver_seq_bj){
            if (!format_out_flag){
               printf("\nBLOCK JACOBI\n");
            }
            SEQ_BlockJacobiRes(A, D, B, Pard, P, &Vars, Params, &Out_data);
         }
         if (solver_seq_bps){
            if (!format_out_flag){
               printf("\nBLOCK PARALLEL SOUTHWELL\n");
            }
            SEQ_BlockParSouthwellRes(A, 
                                     D, 
                                     B, 
                                     Pard, 
                                     P, 
                                     &Vars, 
                                     Params, 
                                     &Out_data);
         }
         PrintResults(A, Vars, Out_data);
      }
      FreeSolveVars(&Vars);
   }
   return 0;
}
