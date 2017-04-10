#ifndef SEQ_SOLVE_H
#define SEQ_SOLVE_H

#include "Southwell.h"

void SEQ_Southwell(CSC A,
                   SolveVars *Vars,
                   SolveParams Params,
                   SolveData *Out_data);

void SEQ_GaussSeidelRes(CSC A,
                        SolveVars *Vars,
                        SolveParams Params,
                        SolveData *Out_data);

void SEQ_JacobiRes(CSC A,
                   SolveVars *Vars,
                   SolveParams Params,
                   SolveData *Out_data);

void SEQ_BlockJacobiRes(CSC A,
                        CSC *D,
                        CSC *B,
                        PardisoInfo *Pard,
                        OrderInfo P,
                        SolveVars *Vars,
                        SolveParams Params,
                        SolveData *Out_data);

void SMEM_JacobiRes(CSC A,
                    SolveVars *Vars,
                    SolveParams Params,
                    SolveData *Out_data);

void SEQ_ParSouthwellRes(CSC A,
                         SolveVars *Vars,
                         SolveParams Params,
                         SolveData *Out_data);

void SMEM_ParSouthwellRes(CSC A,
                          SolveVars *Vars,
                          SolveParams Params,
                          SolveData *Out_data);

void SEQ_BlockParSouthwellRes(CSC A,
                              CSC *D,
                              CSC *B,
                              PardisoInfo *Pard,
                              OrderInfo P,
                              SolveVars *Vars,
                              SolveParams Params,
                              SolveData *Out_data);

void SEQ_MCGS(CSC A,
              OrderInfo P,
              SolveVars *Vars,
              SolveParams Params,
              SolveData *Out_data);

void SMEM_MCGS(CSC A,
               OrderInfo P,
               SolveVars *Vars,
               SolveParams Params,
               SolveData *Out_data);


#endif
