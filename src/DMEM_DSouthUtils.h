#ifndef DMEM_DSOUTH_UTILS_H
#define DMEM_DSOUTH_UTILS_H

#include "Southwell.h"

double DMEM_DSOUTH_RecvDiffNorm2(DMEM_ProcInfo *Proc, int q, int k);

void DMEM_DSOUTH_UpdateResEstim_Norm2(DMEM_MatrixInfo Mat,
                                      DMEM_ProcInfo *Proc,
                                      SolveData *Out_data);

void DMEM_RMA_DSOUTH_GatherBoundRes(DMEM_MatrixInfo Mat,
                                    DMEM_ProcInfo *Proc,
                                    DMEM_SolveVars Vars,
                                    SolveData *Out_data,
                                    int q);

void DMEM_RMA_DSOUTH_RecvResBound(DMEM_MatrixInfo Mat,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveParams Params,
                                  int q,
                                  int k);

void DMEM_RMA_DSOUTH_RecvAllResBound(DMEM_MatrixInfo Mat,
                                     DMEM_ProcInfo *Proc,
                                     DMEM_SolveVars *Vars,
                                     SolveParams Params);

void DMEM_RMA_DSOUTH_RecvAllBound(DMEM_MatrixInfo Mat,
                                  DMEM_ProcInfo *Proc,
                                  DMEM_SolveVars *Vars,
                                  SolveParams Params);

void DMEM_RMA_DSOUTH_AccumExplicitRes(DMEM_ProcInfo *Proc,
                                      SolveData *Out_data,
                                      int q);

void DMEM_POS_DSOUTH_AccumExplicitRes(DMEM_ProcInfo *Proc,
                                      SolveData *Out_data,
                                      int q);

void DMEM_POS_DSOUTH_AccumulateExplicitRes(DMEM_ProcInfo *Proc,
                                           int q);

void DMEM_DSOUTH_ExplicitResUpdate(DMEM_MatrixInfo Mat,
                                   DMEM_ProcInfo *Proc,
                                   SolveData *Out_data,
                                   DMEM_SolveVars *Vars,
                                   SolveParams Params,
                                   VSLStreamStatePtr stream);

void DMEM_DSOUTH_SOS_ExplicitResUpdate(DMEM_MatrixInfo Mat,
                                       DMEM_ProcInfo *Proc,
                                       SolveData *Out_data,
                                       DMEM_SolveVars *Vars,
                                       SolveParams Params,
                                       VSLStreamStatePtr stream);

void DMEM_DSOUTH_OffDiagSweep(DMEM_MatrixInfo Mat,
                              DMEM_ProcInfo *Proc,
                              SolveData *Out_data,
                              SolveVars *Vars);

#endif
