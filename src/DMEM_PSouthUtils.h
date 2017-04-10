#ifndef DMEM_PSOUTH_UTILS_H
#define DMEM_PSOUTH_UTILS_H

#include "Southwell.h"

void DMEM_RMA_PSOUTH_RecvResNorm(DMEM_ProcInfo *Proc);

void DMEM_RMA_PSOUTH_RecvBound(DMEM_MatrixInfo Mat,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars);

void DMEM_SOS_PSOUTH_PutResNormAll(DMEM_ProcInfo *Proc,
                                   SolveData *Out_data);

#endif
