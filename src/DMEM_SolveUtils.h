#ifndef DMEM_SOLVE_UTILS_H
#define DMEM_SOLVE_UTILS_H

#include "Southwell.h"

int DMEM_COLLECT_CheckConverge(int n,
                               DMEM_ProcInfo *Proc,
                               DMEM_SolveVars *Vars,
                               SolveData *Out_data,
                               SolveParams Params);

int DMEM_SOS_CheckConverge(double norm,
                           DMEM_ProcInfo *Proc,
                           SolveData *Out_data,
                           SolveParams Params);

void DMEM_DiagSweep_GS(DMEM_MatrixInfo Mat,
                       DMEM_ProcInfo Proc,
                       SolveData *Out_data,
                       SolveVars *Vars);

void DMEM_OffDiagSweep(DMEM_MatrixInfo Mat,
                       DMEM_ProcInfo *Proc,
                       SolveData *Out_data,
                       SolveVars *Vars);

void DMEM_DiagSweep(DMEM_MatrixInfo Mat,
                    DMEM_ProcInfo *Proc,
                    SolveData *Out_data,
                    DMEM_SolveVars *Vars);

void DMEM_RMA_MCGS_RecvBound(DMEM_MatrixInfo Mat,
                             DMEM_ProcInfo *Proc,
                             double *r,
                             int disp_loc);

void DMEM_RMA_RecvBound(DMEM_MatrixInfo Mat,
                        DMEM_ProcInfo *Proc,
                        double *r);

void DMEM_POS_AccumBound(DMEM_ProcInfo *Proc,
                         SolveData *Out_data,
                         int q);

void DMEM_POS_AccumBoundAll(DMEM_ProcInfo *Proc,
                            SolveData *Out_data);

int DMEM_POS_CheckConverge(double norm,
                           SolveData Out_data,
                           SolveParams Params,
                           DMEM_ProcInfo *Proc);

int DMEM_POS_CheckConverge_array(double norm,
                                 SolveData Out_data,
                                 SolveParams Params,
                                 DMEM_ProcInfo *Proc);

void DMEM_SOS_AccumBoundAll(DMEM_ProcInfo *Proc,
                            SolveData *Out_data);

void DMEM_SOS_Recv(DMEM_ProcInfo *Proc,
                   SolveData *Out_data);

void DMEM_RMA_AccumBound(DMEM_ProcInfo *Proc,
                         int q);
#endif
