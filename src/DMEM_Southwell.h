#ifndef DMEM_SOUTHWELL_H
#define DMEM_SOUTHWELL_H

#include <mpi.h>
#include "Southwell.h"

#define DMEM_SOLVER_J 0
#define DMEM_SOLVER_MCGS 1
#define DMEM_SOLVER_PS 2
#define DMEM_SOLVER_DS 3

#define DMEM_PT2PT 0
#define DMEM_POS 1
#define DMEM_SOS 2

#define DMEM_LOC_SOLVER_DIRECT 0
#define DMEM_LOC_SOLVER_SEQGS 1

typedef struct{
   double **send;
   double **send_prev;
   int *send_disp;
   int *send_count;
   double *recv;
   double *recv_prev;
   int *recv_disp;
   int recv_count;
   MPI_Win win;
}DMEM_RmaBoundInfo;

typedef struct{
   int *send_disp;
   double *recv;
   double *recv_prev;
   MPI_Win win;
}DMEM_RmaResNormInfo;

typedef struct{
   double *norm;
   double *norm_prev;
   double *norm_send_prev;
   double *norm_estim_squ;
   double *my_norm_estim_squ;
   double *my_norm_estim_squ_prev;
   double *neighb_norm_estim;
   int *neighb_explicit_flag;
   int explicit_flag;
   DMEM_RmaResNormInfo Rma;
}DMEM_ResNormInfo;

typedef struct{
   int *num_send_points;
   int *num_recv_points;
   int **send_points;
   int **recv_points;
   int **points;
   int *range;
   int *row_min;
   int *row_max;
   int *row_range_minmax;
   int *row_range_disp;
   double **res;
   double res_norm;
  // DMEM_ResNormInfo Res;
   DMEM_RmaBoundInfo Rma;
}DMEM_BoundInfo;

typedef struct{
   int phase;
   int neighb_flag;
   int my_flag;
   int *recv;
   int win_size;
   MPI_Win win;
}DMEM_ConvergeInfo;

typedef struct{
   /* Number of process neighbors. */
   int size;
   /* Ranks of those neighbors. */
   int *ranks;
   MPI_Group mpi_group;
}DMEM_NeighbInfo;

typedef struct{
   int rank;
   int world_size;
   double res_norm;
   double res_norm_init;
   double res_norm_prev;
   /* Information about MPI process neighbors. */
   DMEM_NeighbInfo Neighb;
   /* Information about boundary points. */
   DMEM_BoundInfo Bound;
   /* Neighborhood residual information. */
   DMEM_ResNormInfo Res;
   DMEM_ConvergeInfo Conv;
}DMEM_ProcInfo;

typedef struct{
   int n_glob;
   int nnz_glob;
   int n;
   int nnz;
   int disp;
   CSC A;
   CSC D;
   CSC *B;
   OrderInfo P;
   PardisoInfo Pard;
}DMEM_MatrixInfo;

typedef struct{
   double *recv;
   int *recv_disp;
   int recv_count;
   int recv_disp_size;
   MPI_Win win;
}DMEM_RmaColorInfo;

typedef struct{
   int n_glob;
   int nnz_glob;
   int n;
   int *disp_loc;
   MPI_Group mpi_group;
   OrderInfo  P_glob;
   OrderInfo *P_loc;
   DMEM_RmaColorInfo *Rma;
   CSC *D;
   CSC ***B;
   DMEM_MatrixInfo *Mat;
   DMEM_NeighbInfo **NeighbSend;
   DMEM_NeighbInfo **NeighbRecv;
   DMEM_BoundInfo **Bound;
}DMEM_MulticolorInfo;

typedef struct{
   SolveVars Loc;
   SolveVars Glob;
}DMEM_SolveVars;

extern int solver_flag;
extern int loc_solver_flag;
extern int implement_flag;
extern int ds_res_estim_flag;
extern int num_samps;

#endif
