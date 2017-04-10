#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include "Southwell.h"

void ReadText_fscanf_metis(FILE *mat_file, MetisGraph *G);

void ReadBinary_fread_metis(FILE *mat_file, MetisGraph *G, Triplet *T);

void Write_csc(CSC A, int base);

void Write_triplet(Triplet T, int base);

void Write_metis(FILE *out_file, MetisGraph G, int base);

void Print_csc(CSC A);

void Print_metis(MetisGraph G);

void ReadText_fscanf_array(char *file_name, double **v, int n);

void ReadVector(char *file_name, double **v, int n);

void WriteBlocks_csc(char *buffer,
                     CSC *D,
                     CSC *B,
                     OrderInfo P,
                     int base);

#endif
