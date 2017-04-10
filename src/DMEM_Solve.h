#ifndef DMEM_SOLVE_H
#define DMEM_SOLVE_H

#include "Southwell.h"

void DMEM_POS_AsyncDistrSouthwell(DMEM_MatrixInfo Mat,
                                  SolveParams Params,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveData *Out_data);

void DMEM_SOS_SyncDistrSouthwell(DMEM_MatrixInfo Mat,
                                 SolveParams Params,
                                 DMEM_ProcInfo *Proc,
                                 DMEM_SolveVars *Vars,
                                 SolveData *Out_data);

void DMEM_POS_AsyncParSouthwell(DMEM_MatrixInfo Mat,
                                SolveParams Params,
                                DMEM_ProcInfo *Proc,
                                DMEM_SolveVars *Vars,
                                SolveData *Out_data);

void DMEM_SOS_SyncParSouthwell(DMEM_MatrixInfo Mat,
                               SolveParams Params,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data);

void DMEM_POS_AsyncBlockJacobi(DMEM_MatrixInfo Mat,
                               SolveParams Params,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data);

void DMEM_SOS_SyncBlockJacobi(DMEM_MatrixInfo Mat,
                              SolveParams Params,
                              DMEM_ProcInfo *Proc,
                              DMEM_SolveVars *Vars,
                              SolveData *Out_data);

void DMEM_SOS_MulticolorGS(DMEM_MulticolorInfo *MC_info,
                           DMEM_ProcInfo Proc,
                           SolveParams Params,
                           DMEM_SolveVars *Vars,
                           SolveData *Out_data);

#endif 
