#ifndef MISC_H
#define MISC_H

#include "Southwell.h"

int InsertInt(int insert, int **x, int n, int max);

double SumAbsDouble(double *x, int n);

double Norm1(double *x, int n);

double SumSquaredDouble(double *x, int n);

double Norm2(double *x, int n);

double SumSquaredAtomicDouble(double *x, int n);

double SumAtomicInt(int *x, int n);

double Norm2Atomic(double *x, int n);

double DiffNorm2Atomic(double *x, double *y, int n);

int Round(double x);

double SumDouble(double *x, int n);

int SumInt(int *x, int n);

unsigned long long SumUL(unsigned long long *x, int n);

double MinMaxScale(double *x, int n, int ind);

double MeanProb(double *x, int n, int ind);

double MaxAbsDouble(double *x, int n);

int IsMaxAbsDouble(double *x, int n, int ind);

int IsZeroMaxAbsDouble(double *x, int n);

int IndMaxAbsDouble(double *x, int n);

double MaxDouble(double *x, int n);

int MaxInt(int *x, int n);

double MinDouble(double *x, int n);

int MinInt(int *x, int n);

double MeanInt(int *x, int n);

double MeanDouble(double *x, int n);

unsigned long long MaxUL(unsigned long long *x, int n);

unsigned long long MinUL(unsigned long long *x, int n);

void SubtractVecDouble(double *x, double *y, int n, double *diff);

void Quicksort_int(int *x, int left, int right);

void QuicksortPair_int_int(int *x, int *y, int left, int right);

void QuicksortPair_int_dbl(int *x, double *y, int left, int right);

void RandDouble(double *v,
                int n,
                double low,
                double high);

void MatVecProd(CSC A, double *x, double *y);

void Residual(CSC A, double *x, double *b, double *r);

void PrintResults(CSC A,
                  SolveVars Vars,
                  SolveData Out_data);

void SetupOutData(SolveData *Out_data, int n);

void ZeroOutData(SolveData *Out_data, int i);

void FreeCSC(CSC *A);

void FreeMetis(MetisGraph *G);

void FreeOrdering(OrderInfo *P);

void FreeSolveVars(SolveVars *Vars);

void DegreeCSC(CSC A, int *degree);

bool AlmostEqual2sComplement_FLT(float A, float B, int maxUlps);

bool AlmostEqual2sComplement_DBL(double A, double B, int maxUlps);

#endif
