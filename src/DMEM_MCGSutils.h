#ifndef DMEM_MCGS_UTILS_H
#define DMEM_MCGS_UTILS_H

#include "Southwell.h"

void DMEM_RMA_MCGS_RecvBound(DMEM_MulticolorInfo *MC_info,
                             DMEM_ProcInfo Proc,
                             double *r,
                             int c);

void DMEM_MCGS_DiagSweep_GS(DMEM_MulticolorInfo MC_info,
                            DMEM_ProcInfo Proc,
                            SolveData *Out_data,
                            SolveVars *Vars,
                            int c);

void DMEM_MCGS_OffDiagSweep_GS(DMEM_MulticolorInfo MC_info,
                               DMEM_ProcInfo Proc,
                               SolveData *Out_data,
                               SolveVars *Vars,
                               int c);

void DMEM_SOS_MCGS_AccumBoundAll(DMEM_MulticolorInfo *MC_info,
                                 DMEM_ProcInfo Proc,
                                 SolveData *Out_data,
                                 int c);

void DMEM_MC_ZeroData(DMEM_MulticolorInfo *MC_info);

#endif
