#include <omp.h>

#include "Southwell.h"
#include "SEQ_Solve.h"
#include "SMEM_Solve.h"
#include "SMEM_Setup.h"
#include "SMEM_Misc.h"
#include "Misc.h"
#include "FileUtils.h"

int mat_file_flag = 0;
int threads = 1;
int format_out_flag = 0;
int color_flag = 0;
int reorder_flag = 1;

int main(int argc, char *argv[])
{
   unsigned long long comm;
   int nnz, n;
   double low = -1, high = 1, rand_range = high-low;
   double rand_num;
   double r_norm;
   double x_low=-1, x_high=1, b_low=-1, b_high=1;
   int num_p = 1, rank_p = 0;

   
   int delay_relax=0;
   int x_file_flag=0, b_file_flag=0;
   int x_zeros_flag=0, x_ones_flag=0;
   int b_zeros_flag=0, b_ones_flag=0;
   int solver_seq_s=0, solver_seq_gs=0;
   int solver_omp_sj=0, solver_omp_sps=0;
   int solver_omp_aj=0, solver_omp_aps=0;
   int solver_omp_mcgs=0;
   unsigned long long aps_pos_comm, abj_pos_comm; 
   unsigned long long agslr_pos_comm, aslr_pos_comm, sj_pos_comm;
   unsigned long long omp_aj_comm, omp_aps_comm, omp_sj_comm, omp_sps_comm;
   int comm_out_flag = 0, time_out_flag = 0;

   double ompwt_start, ompwt_end;
   
   SolveParams InParams; 

   InParams.tol = 1e-5;
   InParams.sweep_max = (int)1e9;
   int arg_iter = 0;
   int m1 = 10, m2 = 10;
   
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
         else if (strcmp(argv[arg_iter], "omp_mcgs") == 0){
            solver_omp_mcgs = 1;
            color_flag = 1;
         }
      }
      else if (strcmp(argv[arg_iter], "-threads") == 0){
         arg_iter++;
         threads = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-tol") == 0){
         arg_iter++;
         InParams.tol = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-sweep_max") == 0){
         arg_iter++;
         InParams.sweep_max = atof(argv[arg_iter]);
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
      else if (strcmp(argv[arg_iter], "-laplace2D_FD5pt_size") == 0){
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
   
   FILE *mat_file = fopen(mat_str, "rb");
   FILE *metis_out_file;

   MetisGraph G;
   SMEM_MatrixInfo Mat;
   SMEM_ThreadInfo Threads;
   SMEM_Setup(mat_file, &Mat, &Threads, m1, m2);
   return 0;
   n = Mat.A.n;
   nnz = Mat.A.nnz;

   SolveData out_seq_s, out_seq_gs, out_omp_aj, out_omp_aps, out_omp_sj, out_omp_sps, out_omp_mcgs;

   if (solver_seq_s){
      out_seq_s.t = (double *)malloc(threads * sizeof(double));
      out_seq_s.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_s.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_s.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
   }

   if (solver_seq_gs){
      out_seq_gs.t = (double *)malloc(threads * sizeof(double));
      out_seq_gs.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_gs.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_seq_gs.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
   }
 
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

   if (solver_omp_mcgs){
      out_omp_mcgs.t = (double *)malloc(threads * sizeof(double));
      out_omp_mcgs.sweep = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mcgs.relax = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mcgs.comm = (unsigned long long *)calloc(threads, sizeof(unsigned long long));
      out_omp_mcgs.relax_hist = (unsigned long long *)calloc(n, sizeof(unsigned long long));
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
   if (solver_seq_s){
      if ((rank_p == 0) && !format_out_flag)
         printf("SEQUENTIAL SOUTHWELL\n");
      ompwt_start = omp_get_wtime();
      SEQ_Southwell(Mat.A, x, b, &u, InParams, &out_seq_s);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_seq_s.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_seq_s.sweep, threads),
                SumUL(out_seq_s.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_seq_s.sweep, threads));
      }
   }
   if (solver_seq_gs){
      if ((rank_p == 0) && !format_out_flag)
         printf("SEQUENTIAL GAUSS-SEIDEL\n");
      ompwt_start = omp_get_wtime();
      SEQ_GaussSeidelRes(Mat.A, x, b, &u, InParams, &out_seq_gs);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_seq_gs.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_seq_gs.sweep, threads),
                SumUL(out_seq_gs.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_seq_gs.sweep, threads));
      }
   }
   if (solver_omp_aps){
      if ((rank_p == 0) && !format_out_flag)
         printf("OpenMP ASYNCHRONOUS PARALLEL SOUTHWELL\n");
      ompwt_start = omp_get_wtime();
      OMP_AsyncParSouthwell(Mat, x, b, &u, InParams, &out_omp_aps);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_omp_aps.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_omp_aps.sweep, threads),
                SumUL(out_omp_aps.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_omp_aps.sweep, threads));
      }
   }
   if (solver_omp_sps){
      if ((rank_p == 0) && !format_out_flag)
         printf("OpenMP SYNCHRONOUS PARALLEL SOUTHWELL\n");
      ompwt_start = omp_get_wtime();
      OMP_SyncParSouthwell(Mat, x, b, &u, InParams, &out_omp_sps);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu %llu\n",
                MaxDouble(out_omp_sps.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_omp_sps.sweep, threads),
                SumUL(out_omp_sps.relax, threads),
                SumUL(out_omp_sps.relax_hist, n));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_omp_sps.sweep, threads));
      }
   }
   if (solver_omp_sj){
      if ((rank_p == 0) && !format_out_flag)
         printf("OpenMP SYNCHRONOUS JACOBI\n");
      ompwt_start = omp_get_wtime();
      OMP_SyncJacobi(Mat, x, b, &u, InParams, &out_omp_sj);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_omp_sj.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_omp_sj.sweep, threads),
                SumUL(out_omp_sj.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_omp_sj.sweep, threads));
      }
   }

   if (solver_omp_aj){
      if ((rank_p == 0) && !format_out_flag)
         printf("OpenMP ASYNCHRONOUS JACOBI\n");
      ompwt_start = omp_get_wtime();
      OMP_AsyncJacobi(Mat, x, b, &u, InParams, &out_omp_aj);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_omp_aj.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_omp_aj.sweep, threads),
                SumUL(out_omp_aj.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_omp_aj.sweep, threads));
      }
   }
   
   if (solver_omp_mcgs){
      if ((rank_p == 0) && !format_out_flag)
         printf("OpenMP MULTICOLOR GAUSS-SEIDEL\n");
      ompwt_start = omp_get_wtime();
      OMP_MCGS(Mat, x, b, &u, InParams, &out_omp_mcgs);
      ompwt_end = omp_get_wtime();
      OMP_Residual(Mat, u, b, &r);
      if (format_out_flag){
         printf("%.10e %.10e %llu %llu\n",
                MaxDouble(out_omp_mcgs.t, threads),
                OMP_Norm2(r, n),
                MaxUL(out_omp_mcgs.sweep, threads),
                SumUL(out_omp_mcgs.relax, threads));
      }
      else {
         printf("time = %f, res norm = %e, sweeps = %lu\n",
                ompwt_end-ompwt_start, Norm2(r, n), MaxUL(out_omp_mcgs.sweep, threads));
      }
   }

   if (!format_out_flag){  
      if (solver_omp_aps && solver_omp_aj)
         if (rank_p == 0)
            printf("\nasync comm ratio = %f\n",
                   (double)SumUL(out_omp_aps.relax,threads)/(double)SumUL(out_omp_aj.relax,threads));
      
      if (solver_omp_sj && solver_omp_sps)
         if (rank_p == 0)
            printf("\nsync comm ratio = %f\n",
                   (double)SumUL(out_omp_sps.relax, threads)/(double)SumUL(out_omp_sj.relax,threads)); 

      if (solver_omp_mcgs && solver_omp_sps)
         if (rank_p == 0)
            printf("\ncomm ratio = %f\n",
                   (double)SumUL(out_omp_sps.relax, threads)/(double)SumUL(out_omp_mcgs.relax,threads));
   }

   
//   free(ones);
//   CSC_Free(&Mat);
//   free(x);
//   free(r);
//   free(b);
//   free(u);
   return 0;
}
