#ifndef SEQ_SETUP_H
#define SEQ_SETUP_H

#include "Southwell.h"

void SEQ_Setup(FILE *in_file,
               CSC *A,
               CSC **B,
               CSC **D,
               OrderInfo *P,
               PardisoInfo **Pard,
               int nx,
               int ny);

void SEQ_SetupNeighb(CSC A, OrderInfo P, SEQ_NeighbInfo *N, int tid);

void SEQ_SetupNeighbBlocks(OrderInfo P,
                           CSC *B,
                           SEQ_BlockNeighbInfo *Neighb);

#endif
