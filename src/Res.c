#define _BSD_SOURCE
#define _XOPEN_SOURCE

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

#include "Misc.h"
#include "IterMethods.h"
#include "List.h"
#include "CSC.h"
#include "MetisGraph.h"
#include "OMP_IterMethods.h"
#include "OMP_Misc.h"
#include "csparse.h"

int mat_file_flag = 0;
int threads = 1;
int format_out_flag = 0;
int color_flag = 0;
int reorder_flag = 1;

int main(int argc, char *argv[])
{
   unsigned long long comm;
   int nnz, np, n;
   double low = -1, high = 1, rand_range = high-low;
   double rand_num;
   double x_low=-1, x_high=1, b_low=-1, b_high=1;
   int num_p = 1, rank_p = 0;

   int delay_relax=0;
   int x_file_flag=0, b_file_flag=0;
   int x_zeros_flag=0, x_ones_flag=0;
   int b_zeros_flag=0, b_ones_flag=0;
   int solver_seq_gs=0, solver_seq_s=0;
   int solver_sj=0;
   int solver_omp_sj=0, solver_omp_sps=0;
   int solver_omp_aj=0, solver_omp_aps=0;
   int solver_omp_mgs=0;
   unsigned long long aps_pos_comm, abj_pos_comm; 
   unsigned long long agslr_pos_comm, aslr_pos_comm, sj_pos_comm;
   unsigned long long omp_aj_comm, omp_aps_comm, omp_sj_comm, omp_sps_comm;
   int comm_out_flag = 0, time_out_flag = 0;

   double mpiwt_start, mpiwt_end;
   MetisGraph G;
   CSC A;

   double tol = 1e-5;
   int arg_iter = 0;
   int m1 = 10, m2 = 10;
   double t_low = -4, t_high = -4, t_iter = -4;
   unsigned long long max_sweep = (int)1e9;
   unsigned long long max_comm = (int)1e9;
   char mat_str[100] = "";
   char x_file_str[100] = "";
   char b_file_str[100] = "";
   while(arg_iter < argc){
      if (strcmp(argv[arg_iter], "-solver") == 0){
         arg_iter++;
         if (strcmp(argv[arg_iter], "seq_s") == 0){
            solver_seq_s = 1;
         }
         else if (strcmp(argv[arg_iter], "seq_gs") == 0){
            solver_seq_gs = 1;
         }
         else if (strcmp(argv[arg_iter], "sj") == 0){
            solver_sj = 1;
         }
         else if (strcmp(argv[arg_iter], "omp_aps") == 0){
            solver_omp_aps = 1;
         }
         else if (strcmp(argv[arg_iter], "omp_aj") == 0){
            solver_omp_aj = 1;
         }
         else if (strcmp(argv[arg_iter], "omp_sps") == 0){
            solver_omp_sps = 1;
         }
         else if (strcmp(argv[arg_iter], "omp_sj") == 0){
            solver_omp_sj = 1;
         }
         else if (strcmp(argv[arg_iter], "omp_mgs") == 0){
            solver_omp_mgs = 1;
            color_flag = 1;
         }
      }
      else if (strcmp(argv[arg_iter], "-threads") == 0){
         arg_iter++;
         threads = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-t_low") == 0){
         arg_iter++;
         t_low = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-t_high") == 0){
         arg_iter++;
         t_high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-t_iter") == 0){
         arg_iter++;
         t_iter = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-tol") == 0){
         arg_iter++;
         tol = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-max_sweep") == 0){
         arg_iter++;
         max_sweep = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-max_comm") == 0){
         arg_iter++;
         max_comm = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-delay") == 0){
         arg_iter++;
         delay_relax = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-mat_file") == 0){
         arg_iter++;
         strcpy(mat_str, argv[arg_iter]);
         mat_file_flag = 1;
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
      else if (strcmp(argv[arg_iter], "-b_rand_range") == 0){
          arg_iter++;
          b_low = atof(argv[arg_iter++]);
          b_high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-fd_2d_size") == 0){
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

   if ((solver_seq_gs || solver_seq_s) && (threads > 1)){
      printf("Error: cannot more than one thread with sequential Southwell or Gauss-Seidel\n");
      exit(1);
   }
   
   FILE *mat_file = fopen(mat_str, "r");
   FILE *metis_out_file;

   Neighbors neighb;
   CSC_Init(mat_file, &A, &neighb, m1, m2);
   n = A.n;
   np = A.np[rank_p];
   nnz = A.nnz;

   OutData out_omp_aj, out_omp_aps, out_omp_sj, out_omp_sps, out_omp_mgs;
   OutData out_seq_s, out_seq_gs;
 
   if (solver_omp_aj){
      out_omp_aj.t = (double *)malloc(threads * sizeof(double));
      out_omp_aj.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_aj.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_aj.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
   }
   
   if (solver_omp_aps){
      out_omp_aps.t = (double *)malloc(threads * sizeof(double));
      out_omp_aps.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_aps.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_aps.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_aps.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
   }

   if (solver_omp_sj){
      out_omp_sj.t = (double *)malloc(threads * sizeof(double));
      out_omp_sj.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_sj.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_sj.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
   }

   if (solver_omp_sps){
      out_omp_sps.t = (double *)malloc(threads * sizeof(double));
      out_omp_sps.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_sps.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_sps.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_sps.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
   }

   if (solver_omp_mgs){
      out_omp_mgs.t = (double *)malloc(threads * sizeof(double));
      out_omp_mgs.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mgs.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mgs.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mgs.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
   }

   if (solver_seq_gs){
      out_seq_gs.t = (double *)malloc(threads * sizeof(double));
      out_seq_gs.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_gs.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_gs.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_gs.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
   }

   if (solver_seq_s){
      out_seq_s.t = (double *)malloc(threads * sizeof(double));
      out_seq_s.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_s.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_s.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_s.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
   }

   
   double *r = (double *)calloc(n, sizeof(double));
   double *x = (double *)calloc(n, sizeof(double));
   double *b = (double *)calloc(n, sizeof(double));
   double *u = (double *)calloc(n, sizeof(double));
   double *ones = (double *)calloc(n, sizeof(double));
   for (int i = 0; i < n; i++) ones[i] = 1.0;
   srand((unsigned)time(NULL));

   if (x_file_flag){
      ReadVector(x_file_str, &x, n);
   }
   else if (x_zeros_flag){
      for (int i = 0; i < n; i++) x[i] = 0;
   }
   else if (x_ones_flag){
      for (int i = 0; i < n; i++) x[i] = 1.0;
   }
   else {
      for (int i = 0; i < n; i++) x[i] = RandDouble(x_low, x_high);
   }

   if (b_file_flag){
      ReadVector(b_file_str, &b, n);
   }
   else if (b_zeros_flag){
      for (int i = 0; i < n; i++) b[i] = 0;
   }
   else if (b_ones_flag){
      for (int i = 0; i < n; i++) b[i] = 1.0;
   }
   else {
      for (int i = 0; i < n; i++) b[i] = RandDouble(b_low, b_high);
   }
   for (int i = 0; i < n; i++){ 
      u[i] = x[i];
      ones[i] = 1;
   }


   /*************************
    * OpenMP
    * ***********************/ 
   double t = t_low;
   double t_obsvd = 0, t_obsvd_prev = -1;
   unsigned long long sweep, sweep_prev = 0;
   double r_norm, r_norm_prev = DBL_MAX;
   while (t >= t_high){
      tol = pow(10, t);
      if (solver_omp_aps){
         OMP_AsyncParSouthwell(A, x, b, &u, tol, max_sweep, &out_omp_aps);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         t_obsvd = MeanDouble(out_omp_aps.t, threads);
         sweep = MinUL(out_omp_aps.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MeanDouble(out_omp_aps.t, threads),
                   r_norm,
                   MinUL(out_omp_aps.sweep, threads),
                   SumUL(out_omp_aps.relax, threads));
         }
      }
      if (solver_omp_sps){
         OMP_SyncParSouthwell(A, x, b, &u, tol, max_sweep, &out_omp_sps);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         sweep = MinUL(out_omp_sps.sweep, threads);
         t_obsvd = MeanDouble(out_omp_sps.t, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu %llu\n",
                   MeanDouble(out_omp_sps.t, threads),
                   r_norm,
                   MinUL(out_omp_sps.sweep, threads),
                   SumUL(out_omp_sps.relax, threads),
                   SumUL(out_omp_sps.relax_hist, n));
         }
         
      }
      if (solver_omp_sj){
         OMP_SyncJacobi(A, x, b, &u, tol, max_sweep, &out_omp_sj);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         t_obsvd = MeanDouble(out_omp_sj.t, threads);
         sweep = MinUL(out_omp_sj.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MeanDouble(out_omp_sj.t, threads),
                   r_norm,
                   MinUL(out_omp_sj.sweep, threads),
                   SumUL(out_omp_sj.relax, threads));
         }
      }

      if (solver_omp_aj){
         OMP_AsyncJacobi(A, x, b, &u, tol, max_sweep, &out_omp_aj);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         t_obsvd = MeanDouble(out_omp_aj.t, threads);
         sweep = MinUL(out_omp_aj.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MeanDouble(out_omp_aj.t, threads),
                   r_norm,
                   MinUL(out_omp_aj.sweep, threads),
                   SumUL(out_omp_aj.relax, threads));
         } 
      }
      
      if (solver_omp_mgs){
         OMP_MGS(A, x, b, &u, tol, max_sweep, &out_omp_mgs);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         t_obsvd = MeanDouble(out_omp_mgs.t, threads);
         sweep = MinUL(out_omp_mgs.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MeanDouble(out_omp_mgs.t, threads),
                   r_norm,
                   MinUL(out_omp_mgs.sweep, threads),
                   SumUL(out_omp_mgs.relax, threads));
         }
      }
      if (solver_seq_s){
         Southwell(A, x, b, &u, tol, max_sweep, &out_seq_s);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         sweep = MinUL(out_seq_s.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MaxDouble(out_seq_s.t, threads),
                   r_norm,
                   MinUL(out_seq_s.sweep, threads),
                   SumUL(out_seq_s.relax, threads));
         }
      }
      if (solver_seq_gs){
         GaussSeidelRes(A, x, b, &u, tol, max_sweep, &out_seq_gs);
         OMP_Residual(A, u, b, &r);
         r_norm = OMP_Norm2(r, n);
         sweep = MinUL(out_seq_gs.sweep, threads);
         if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
            printf("%.10e %.10e %llu %llu\n",
                   MaxDouble(out_seq_gs.t, threads),
                   r_norm,
                   MinUL(out_seq_gs.sweep, threads),
                   SumUL(out_seq_gs.relax, threads));
         }
      }
      if ((sweep > sweep_prev) && (r_norm < r_norm_prev)){
         t_obsvd_prev = t_obsvd;
         sweep_prev = sweep;
         r_norm_prev = r_norm;
      }
      t += t_iter;
   }
   
//   free(ones);
//   CSC_Free(&A);
//   free(x);
//   free(r);
//   free(b);
//   free(u);
   return 0;
}
