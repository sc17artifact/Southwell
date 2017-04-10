#include "Southwell.h"
#include "Misc.h"

int InsertInt(int insert, int **x, int n, int max)
{
   for (int i = 0; i < n; i++){
      if ((*x)[i] == insert) return 0;
      else if ((*x)[i] >= max){
         (*x)[i] = insert;
         return 1;
      }
   }
   return 0;
}

/*
 * absolute value sum of a double array.
 */
double SumAbsDouble(double *x, int n)
{
   double sum = 0;
   for (int i = 0; i < n; i++) sum += fabs(x[i]);
   return sum;
}

/*
 * compute the 1-norm of an array.
 */
double Norm1(double *x, int n)
{
   double sum = 0;
   for (int i = 0; i < n; i++){
      sum += fabs(x[i]);
   }
   return sum;
}

double SumSquaredDouble(double *x, int n)
{
   double sum = 0;
   for (int i = 0; i < n; i++){
      sum += pow(fabs(x[i]), 2);
   }
   return sum;
}

//double SumAtomicInt(int *x, int n)
//{
//   int sum = 0, elem;
//   for (int i = 0; i < n; i++){
//      #pragma omp atomic read
//      elem = x[i];
//      sum += elem;
//   }
//   return sum;
//}

//double SumSquaredAtomicDouble(double *x, int n)
//{
//   double sum = 0, elem;
//   for (int i = 0; i < n; i++){
//      #pragma omp atomic read
//      elem = x[i];
//      sum += pow(fabs(elem), 2);
//   }
//   return sum;
//}

/*
 * compute the 2-norm of an array.
 */
double Norm2(double *x, int n)
{
   return sqrt(SumSquaredDouble(x, n));
}

//double Norm2Atomic(double *x, int n)
//{
//   return sqrt(SumSquaredAtomicDouble(x, n));
//}

//double DiffNorm2Atomic(double *x, double *y, int n)
//{
//   double sum = 0;
//   double x_elem, y_elem;
//   for (int i = 0; i < n; i++){
//      #pragma omp atomic read
//      x_elem = x[i];
//      #pragma omp atomic read
//      y_elem = y[i];
//
//      sum += pow(x_elem - y_elem, 2);
//   }
//   return sqrt(sum);
//}


/*
 * round a double number.
 */
int Round(double x){
    if (x >= 0.0) return (int)(0.5 + x);
    else return (int)(-0.5 + x);
}

/*
 * sum of an array of doubles.
 */
double SumDouble(double *x, int n)
{
   double s = 0;
   for (int i = 0; i < n; i++){
      s += x[i];
   }
   return s;
}

/*
 * sum of an array of integers.
 */
int SumInt(int *x, int n)
{
   int s = 0;
   for (int i = 0; i < n; i++){
      s += x[i];
   }
   return s;
}

unsigned long long SumUL(unsigned long long *x, int n)
{
   unsigned long long s = 0;
   for (int i = 0; i < n; i++){
      s += x[i];
   }
   return s;
}

/*
 * determine the min-max probability of element ind relative to the other 
 * array elements.
 */
double MinMaxScale(double *x, int n, int ind)
{
   double den, max_val = 0, min_val, scaled_val;
   scaled_val = fabs(x[ind]);
 
   if (n > 1){
      for (int i = 0; i < n; i++){
         if (max_val < fabs(x[i])) max_val = fabs(x[i]);
      }
      min_val = max_val;
      for (int i = 0; i < n; i++)
         if (min_val > fabs(x[i])) min_val = fabs(x[i]);

      den = max_val - min_val;
      if (fabs(den) <=  DBL_MIN) return 1; 
      else {
         scaled_val = (scaled_val - min_val) / (max_val - min_val);
         return scaled_val;
      }
   } 
   return 1;
}

/*
 * determine the probability of element ind relative to the other array 
 * elements.
 */
double MeanProb(double *x, int n, int ind)
{
   double sum = 0;
   for (int i = 0; i < n; i++){ 
      sum += fabs(x[i]);
   }
   if (sum <= DBL_MIN) return 0;
   return fabs(x[ind]) / sum;
}

double MaxAbsDouble(double *x, int n)
{
   int i_max = 0;
   double max_val = 0;
   for(int i = 0; i < n; i++){
      if(fabs(x[i]) > fabs(max_val)){
         max_val = x[i];
         i_max = i;
      }
   }
   return fabs(max_val);
}

int IsMaxAbsDouble(double *x, int n, int ind)
{
   int i_max = 0;
   double max_val = 0;
   for(int i = 0; i < n; i++){
      if(fabs(x[i]) > fabs(max_val)){
         max_val = x[i];
         i_max = i;
      }
   }
   if (ind == i_max){
      return 1;
   }
   else{
      return 0;
   } 
}

int IsZeroMaxAbsDouble(double *x, int n)
{
   double val_max = fabs(x[0]);
   double val_check;
   for(int i = 1; i < n; i++){
      val_check = fabs(x[i]);
      if (!AlmostEqual2sComplement_DBL(val_check, val_max, 1)){
         if(val_check > val_max) return 0;
      }
   }
   return 1;
}

/* determine the index of the max absolute value element in an array. */
int IndMaxAbsDouble(double *x, int n)
{
   int i_max = 0;
   double max_val = 0;
   for(int i = 0; i < n; i++){
      if(fabs(x[i]) > fabs(max_val)){
         max_val = x[i];
         i_max = i;
      }
   }
   return i_max;
}

double MaxDouble(double *x, int n)
{
   double max_val = x[0];
   for(int i = 1; i < n; i++){
      if(x[i] > max_val){
         max_val = x[i];
      }
   }
   return max_val;
}

int MaxInt(int *x, int n)
{
   int max_val = 0;
   for(int i = 0; i < n; i++){
      if(x[i] > max_val){
         max_val = x[i];
      }
   }
   return max_val;
}

double MinDouble(double *x, int n)
{
   double min_val = x[0];
   for(int i = 1; i < n; i++){
      if(x[i] < min_val){
         min_val = x[i];
      }
   }
   return min_val;
}

int MinInt(int *x, int n)
{
   int min_val = x[0];
   for(int i = 1; i < n; i++){
      if(x[i] < min_val){
         min_val = x[i];
      }
   }
   return min_val;
}

double MeanInt(int *x, int n)
{
   double sum = 0;
    
   for(int i = 0; i < n; i++){
      sum += (double)x[i];
   }
   return sum/n;
}

double MeanDouble(double *x, int n)
{
   double sum = 0;
   for(int i = 0; i < n; i++){
      sum += (double)x[i];
   }
   return sum/n;
}


unsigned long long MaxUL(unsigned long long *x, int n)
{
   int i_max = 0;
   unsigned long long max_val = 0;
   for(int i = 0; i < n; i++){
      if(x[i] > max_val){
         max_val = x[i];
      }
   }
   return max_val;
}

unsigned long long MinUL(unsigned long long *x, int n)
{
   unsigned long long min_val = x[0];
   for(int i = 1; i < n; i++){
      if(x[i] < min_val){
         min_val = x[i];
      }
   }
   return min_val;
}

void SubtractVecDouble(double *x, double *y, int n, double *diff)
{
   for (int i = 0; i < n; i++) diff[i] = x[i] - y[i];
}

void Quicksort_int(int *x, int left, int right)
{
   int i = left, j = right+1, pivot = x[left], temp;
   if (left < right){
      while(1){
         do{
            ++i;
         }while((x[i] <= pivot) && (i <= right));
         do{
            --j;
         }while(x[j] > pivot);
         if (i >= j) break;
         temp = x[i];
         x[i] = x[j];
         x[j] = temp;
      }
      temp = x[left];
      x[left] = x[j];
      x[j] = temp;
      Quicksort_int(x, left, j-1);
      Quicksort_int(x, j+1, right);
   }
}

void QuicksortPair_int_int(int *x, int *y, int left, int right)
{
   int i = left, j = right+1, pivot = x[left], temp;
   if (left < right){
      while(1){
         do{
            ++i;
         }while((x[i] <= pivot) && (i <= right));
         do{
            --j;
         }while(x[j] > pivot);
         if (i >= j) break;
         temp = x[i];
         x[i] = x[j];
         x[j] = temp;
         temp = y[i];
         y[i] = y[j];
         y[j] = temp;
      }
      temp = x[left];
      x[left] = x[j];
      x[j] = temp;
      temp = y[left];
      y[left] = y[j];
      y[j] = temp;
      QuicksortPair_int_int(x, y, left, j-1);
      QuicksortPair_int_int(x, y, j+1, right);
   }
}

void QuicksortPair_int_dbl(int *x, double *y, int left, int right)
{
   int i = left, j = right+1, pivot = x[left], temp;
   double temp_dbl;
   if (left < right){
      while(1){
         do{
            ++i;
         }while((x[i] <= pivot) && (i <= right));
         do{
            --j;
         }while(x[j] > pivot);
         if (i >= j) break;
         temp = x[i];
         x[i] = x[j];
         x[j] = temp;
         temp_dbl = y[i];
         y[i] = y[j];
         y[j] = temp_dbl;
      }
      temp = x[left];
      x[left] = x[j];
      x[j] = temp;
      temp_dbl = y[left];
      y[left] = y[j];
      y[j] = temp_dbl;
      QuicksortPair_int_dbl(x, y, left, j-1);
      QuicksortPair_int_dbl(x, y, j+1, right);
   }
}

void RandDouble(double *v,
                int n,
                double low,
                double high)
{
   VSLStreamStatePtr stream;
   vslNewStream(&stream, VSL_BRNG_SFMT19937, 0);
   vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, v, low, high);
}

void MatVecProd(CSC A, double *x, double *y)
{
   int row, ind;
   for (int i = 0; i < A.n; i++) y[i] = 0;
   for (int i = 0; i < A.n; i++){
      for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++){
         ind = A.j_ptr[i]+j;
         row = A.i[ind];
         y[row] += A.a[ind] * x[i];
      }
   }
}

double VecVecProd(double *x, double *y, int n)
{
   double sum = 0;
   for (int i = 0; i < n; i++){
      sum += x[i] * y[i];
   }
   return sum;
}

double ErrAnorm(CSC A, double *e)
{
   double *Ae = (double *)calloc(A.n, sizeof(double));
   MatVecProd(A, e, Ae);
   double Aee = VecVecProd(Ae, e, A.n);
   free(Ae);
   return sqrt(Aee);
}


void Residual(CSC A, 
              double *x, 
              double *b,
              double *r)
{
   double *Ax = (double *)calloc(A.n, sizeof(double));
   MatVecProd(A, x, Ax);
   for (int i = 0; i < A.n; i++){
      r[i] = b[i] - Ax[i];
   }
   free(Ax);
}

void DegreeCSC(CSC A, int *degree)
{
   for (int i = 0; i < A.n; i++)
      degree[i] = (A.j_ptr[i+1] - A.j_ptr[i]);
}

void PrintResults(CSC A,
                  SolveVars Vars,
                  SolveData Out_data)
{
   double r_norm;
  // Residual(A, Vars.u, Vars.b, Vars.r);
  // r_norm = Norm2(Vars.r, A.n);
  // for (int i = 0; i < A.n; i++) printf("%f\n", Vars.u[i]);  
   double *ones = (double *)calloc(A.n, sizeof(double));
   for (int i = 0; i < A.n; i++) ones[i] = 1;
   MatVecProd(A, Vars.x, Vars.b);
   double *e = (double *)calloc(A.n, sizeof(double));
   double b_norm = Norm2(Vars.b, A.n);
   for (int i = 0; i < A.n; i++){
      e[i] = Vars.x[i]/b_norm - Vars.u[i];
   }
   r_norm = ErrAnorm(A, e);
   if (format_out_flag){
      printf("%e %e %llu %e %e ",
             r_norm,
             Out_data.wtime_tot[0],
             Out_data.sweep[0],
             Out_data.relax_scaled[0],
             Out_data.work_iter[0]);
      printf("\n");
   }
   else {
      printf("\nnorm = %e, wall-time = %e, sweeps = %llu, (total relax)/n = %e, work = %e\n",
             r_norm,
             Out_data.wtime_tot[0],
             Out_data.sweep[0],
             Out_data.relax_scaled[0],
             Out_data.work_iter[0]);
      printf("\n");
   }
}

void SetupOutData(SolveData *Out_data, int n)
{
   Out_data->wtime_tot = (double *)calloc(n, sizeof(double));
   Out_data->wtime_comp = (double *)calloc(n, sizeof(double));
   Out_data->wtime_comm = (double *)calloc(n, sizeof(double));
   Out_data->wtime_conv = (double *)calloc(n, sizeof(double));
   Out_data->wtime_iter_tot = (double *)calloc(n, sizeof(double));
   Out_data->wtime_iter_comp = (double *)calloc(n, sizeof(double));
   Out_data->wtime_iter_comm = (double *)calloc(n, sizeof(double));
   Out_data->wtime_comm_res = (double *)calloc(n, sizeof(double));
   Out_data->wtime_comm_sweep = (double *)calloc(n, sizeof(double));
   Out_data->wtime_update_bound = (double *)calloc(n, sizeof(double));
   Out_data->wtime_res_estim = (double *)calloc(n, sizeof(double));
   Out_data->comm_scaled = (double *)calloc(n, sizeof(double));
   Out_data->relax_scaled = (double *)calloc(n, sizeof(double));
   Out_data->flop_scaled = (double *)calloc(n, sizeof(double));
   Out_data->sweep =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->relax =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->flop =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->comm =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->comm_res =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->comm_sweep =
      (unsigned long long *)calloc(n, sizeof(unsigned long long));
   Out_data->relax_mask = (int *)calloc(n, sizeof(int));
   Out_data->relax_hist = (int *)calloc(n, sizeof(int));
   Out_data->work_iter = (double *)calloc(n, sizeof(double));
}

void ZeroOutData(SolveData *Out_data, int i)
{
   Out_data->wtime_tot[i] = 0;
   Out_data->wtime_comp[i] = 0;
   Out_data->wtime_comm[i] = 0;
   Out_data->wtime_conv[i] = 0;
   Out_data->wtime_iter_tot[i] = 0;
   Out_data->wtime_iter_comm[i] = 0;
   Out_data->wtime_iter_comp[i] = 0;
   Out_data->wtime_comm_res[i] = 0;
   Out_data->wtime_comm_sweep[i] = 0;
   Out_data->wtime_update_bound[i] = 0;
   Out_data->wtime_res_estim[i] = 0;
   Out_data->sweep[i] = 0;
   Out_data->relax[i] = 0;
   Out_data->relax_scaled[i] = 0;
   Out_data->comm[i] = 0;
   Out_data->comm_res[i] = 0;
   Out_data->comm_sweep[i] = 0;
   Out_data->comm_scaled[i] = 0;
   Out_data->flop[i] = 0;
   Out_data->flop_scaled[i] = 0;
}

void FreeOrdering(OrderInfo *P)
{
   free(P->dispv);
   free(P->disp);
   free(P->part);
}

void FreeMetis(MetisGraph *G)
{
   free(G->adjncy);
   free(G->adjwgt);
   free(G->xadj);
}

void FreeCSC(CSC *A)
{
   free(A->a);
   free(A->i);
   free(A->j_ptr);
   free(A->diag);
}

void FreeSolveVars(SolveVars *Vars)
{
   free(Vars->r);
   free(Vars->u);
   free(Vars->b);
   free(Vars->du);
}

bool AlmostEqual2sComplement_DBL(double A, double B, int maxUlps) 
{ 
   assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024); 
   int64_t aLong = *reinterpret_cast<int64_t*>( &A ); 
   if (aLong < 0){ 
      aLong = 0x8000000000000000 - aLong; 
   }
   int64_t bLong = *reinterpret_cast<int64_t*>( &B ); 
   if (bLong < 0){ 
      bLong = 0x8000000000000000 - bLong; 
   }
   int64_t longDiff = (aLong - bLong) & 0x7FFFFFFFFFFFFFFF; 
   if (longDiff <= maxUlps){ 
      return true; 
   }
   return false; 
}

bool AlmostEqual2sComplement_FLT(float A, float B, int maxUlps)
{
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    int aInt = *(int*)&A;
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;
    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;
    int intDiff = abs(aInt - bInt);
    if (intDiff <= maxUlps)
        return true;
   return false;
}
