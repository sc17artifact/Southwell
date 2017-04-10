#ifndef DMEM_SETUP_H
#define DMEM_SETUP_H

#include "Southwell.h"
#include "DMEM_Southwell.h"

void DMEM_MetisSetup(FILE *in_file,
                     DMEM_MatrixInfo *Mat,
                     DMEM_ProcInfo *Proc,
                     int m,
                     int w);

void DMEM_MulticolorSetup(FILE *in_file,
                          DMEM_MulticolorInfo *MC_info,
                          DMEM_ProcInfo *Proc,
                          int m,
                          int w);

void DMEM_SetupOutData(SolveData *Out_data);

void DMEM_ZeroOutData(SolveData *Out_data);

void DMEM_ZeroProc(DMEM_ProcInfo *Proc);

#endif
