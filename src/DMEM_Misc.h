#ifndef DMEM_MISC_H
#define DMEM_MISC_H

double RmaResNorm2(double *r, int n);

void DMEM_Residual(CSC A, double *x_p, double *b,
                   double **r);

void DMEM_Write_csc(char *buffer,
                    DMEM_MatrixInfo Mat,
                    DMEM_ProcInfo Proc,
                    int base);

void DMEM_WriteBlocks_csc(char *buffer,
                          DMEM_MatrixInfo Mat,
                          DMEM_ProcInfo Proc,
                          int base);

void DMEM_WriteBlocks_MC(char *buffer,
                         DMEM_MulticolorInfo MC_info,
                         DMEM_ProcInfo Proc,
                         int base);
void DMEM_MatVecProd_MC_CscBlocks(DMEM_MulticolorInfo MC_info,
                                  double *x,
                                  double *z,
                                  int n_glob);

void DMEM_DistributeVector_dbl(double *x_glob,
                               double *x,
                               DMEM_ProcInfo Proc,
                               DMEM_MatrixInfo Mat);

void DMEM_MatVecProd_CscBlocks(DMEM_MatrixInfo Mat,
                               DMEM_ProcInfo Proc,
                               double *x,
                               double *z,
                               int n_glob);

void DMEM_Residual_CscBlocks(DMEM_MatrixInfo Mat,
                             DMEM_ProcInfo Proc,
                             double *x,
                             double *b,
                             double *r);

void DMEM_Residual(DMEM_MatrixInfo Mat,
                   DMEM_MulticolorInfo MC_info,
                   DMEM_ProcInfo Proc,
                   DMEM_SolveVars *Vars);

void DMEM_PrintResults(DMEM_MatrixInfo Mat,
                       DMEM_MulticolorInfo MC_info,
                       DMEM_ProcInfo Proc,
                       DMEM_SolveVars *Vars,
                       SolveData Out_data);

void DMEM_RandDouble(DMEM_ProcInfo Proc,
                     double *v,
                     int n,
                     double low,
                     double high);

void DMEM_Free_SolveVars(DMEM_SolveVars *Vars);

void DMEM_Free_MatrixInfo(DMEM_MatrixInfo *Mat,
                         int neighb_size);

void DMEM_Free_ProcInfo(DMEM_ProcInfo *Proc);

double DMEM_WINREAD_DiffNorm2(DMEM_ProcInfo *Proc, int q, int k);

double DMEM_BoundSumSqu(DMEM_ProcInfo *Proc,
                        DMEM_SolveVars Vars);

#endif
