#ifndef SMEM_SOUTHWELL_H
#define SMEM_SOUTHWELL_H

#include "Southwell.h"

#define SMEM_SOLVER_SJ 0
#define SMEM_SOLVER_SPS 1
#define SMEM_SOLVER_AJ 2
#define SMEM_SOLVER_APS 3

typedef struct{
   int *len;
   int **rows;
   double **res;
}SMEM_NeighbInfo;

typedef struct{
   int *interior;
   int *bound;
}SMEM_BoundInfo;

typedef struct{
   SMEM_NeighbInfo neighb;
   SMEM_BoundInfo bound;
}SMEM_ThreadInfo;

typedef struct{
   CSC A;
   OrderInfo P;
}SMEM_MatrixInfo;

typedef struct{
   double *data;
   int *map;
   int *map_loc;
   int low;
   int high;
   int data_len;
   int map_len;
   int iter;
}SMEM_RelaxMap;

extern int threads;

#endif
