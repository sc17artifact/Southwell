#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "Southwell.h"

void ReorderTriplet_csc(Triplet T, CSC *A, OrderInfo *P);

void Reorder(OrderInfo *P,
             Triplet *T,
             CSC *A);

void Metis_to_Triplet(MetisGraph G, Triplet *T);

void Metis_to_CSC(CSC *Mat, MetisGraph G);

void Get_Loc_CSC(CSC A, OrderInfo P, CSC *A_loc, int parti);

void List_to_Block(CSC *A,
                   std::vector<std::list<int>> ind_list,
                   std::vector<std::list<double>> elem_list);

void ScaleDiag_CSC(CSC *A);

#endif
