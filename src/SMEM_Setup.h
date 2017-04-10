#include "Southwell.h"

void SMEM_Setup(
   FILE *in_file,
   SMEM_MatrixInfo *Mat,
   SMEM_ThreadInfo *Threads,
   int m,
   int w);

void OMP_SetupNeighb(SMEM_MatrixInfo A, SMEM_ThreadInfo *T);

void OMP_SetupNeighb_thread(SMEM_MatrixInfo Mat, SMEM_ThreadInfo *T);

void OMP_SetupWrite_scatter(SMEM_MatrixInfo Mat, Write *W);

void OMP_SetupWrite_gather(SMEM_MatrixInfo Mat, Write *W);

