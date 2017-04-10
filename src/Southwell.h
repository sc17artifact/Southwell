#ifndef SOUTHWELL_H
#define SOUTHWELL_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <algorithm>
#include <list>
#include <vector>
#include <mm_malloc.h>
#include <time.h>
#include <functional>
#include <omp.h>
#include <mkl.h>
#include <mkl_vsl.h>

#include "metis.h"

typedef struct{
   double *a;
   int *i;
   int *j_ptr;
   double *diag;
   int n;
   int nnz;
}CSC;

typedef struct{
   double *a;
   MKL_INT *ia;
   MKL_INT *ja;
   MKL_INT n;
   MKL_INT nnz;
}PardisoCSR;

typedef struct{
   double wtime_setup;
   void *pt[64];
   MKL_INT maxfct;
   MKL_INT mnum;
   MKL_INT mtype;
   MKL_INT phase;
   MKL_INT *perm;
   MKL_INT nrhs;
   MKL_INT iparm[64];
   MKL_INT msglvl; 
   MKL_INT error;
   MKL_INT ddum; 
   MKL_INT idum;
   PardisoCSR csr;
}PardisoInfo;

typedef struct{
   int nparts;
   int *dispv;
   int *disp;
   int *part;
   int *perm;
   int *map;
}OrderInfo;

typedef struct{
   int *len;
   double **res;
}SEQ_NeighbInfo;

typedef struct{
   int *size;
   int **blocks;
}SEQ_BlockNeighbInfo;

typedef struct{
   int *i;
   int *j;
   double *a;
   int n;
   int nnz;
}Triplet;

typedef struct{
   int i;
   int j;
   double a;
}Triplet_AOS;

typedef struct{
   idx_t *xadj;
   idx_t *adjncy;
   real_t *adjwgt;
   idx_t n;
   idx_t nnz;
}MetisGraph;

typedef struct{
   double tol;
   int ds_delay;
   unsigned long long sweep_max;
   unsigned long long relax_max;
}SolveParams;

typedef struct{
   double *wtime_tot;
   double *wtime_comp;
   double *wtime_comm;
   double *wtime_conv;
   double *wtime_iter_tot;
   double *wtime_iter_comp;
   double *wtime_iter_comm;
   double *wtime_comm_res;
   double *wtime_comm_sweep;
   double *wtime_update_bound;
   double *wtime_res_estim;
   double *relax_scaled;
   double *comm_scaled;
   double *flop_scaled;
   unsigned long long *sweep;
   unsigned long long *flop;
   unsigned long long *relax;
   unsigned long long *comm;
   unsigned long long *comm_res;
   unsigned long long *comm_sweep;
   int *relax_hist;
   int *relax_mask;
   double *work_iter;
}SolveData;

typedef struct{
   int ds_delay_count;
   double *x;
   double *b;
   double *r;
   double *u;
   double *du;
}SolveVars;

extern int mat_file_flag;
extern int color_flag;
extern int format_out_flag;
extern int symmetric_flag;

#endif
