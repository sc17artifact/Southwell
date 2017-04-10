#include "Southwell.h"
#include "Misc.h"
#include "SEQ_Setup.h"

extern int threads;
extern int loc_direct_flag;
extern int ps_num_relax;

void BlockSweepGS(CSC A,
                  OrderInfo P,
                  SolveVars *Vars,
                  SolveData *Out_data,
                  double **r_loc,
                  int p);

void BlockSweepDirect(CSC D,
                      CSC B,
                      PardisoInfo Pard,
                      OrderInfo P,
                      SolveVars *Vars,
                      SolveData *Out_data,
                      double **r_loc,
                      int p);

void BlockOffDiagSweepDirect(CSC *B,
                             OrderInfo P,
                             SolveVars *Vars,
                             SolveData *Out_data,
                             double **r_loc);

double OMP_Norm2(double *x,
                 int n)
{
   double sum = 0;
   #pragma omp parallel for num_threads(threads) reduction(+:sum)
   for (int i = 0; i < n; i++){
      sum += pow(fabs(x[i]), 2);
   }
   #pragma omp barrier
   return sqrt(sum);
}

void SEQ_Southwell(CSC A, 
                   SolveVars *Vars,
                   SolveParams Params,
                   SolveData *Out_data)
{
   int i_max, ind, row;
   double r_norm, du, u_max, elem, comm = 0, relax = 0;
   double start = omp_get_wtime();
   while((Norm2(Vars->r, A.n) > Params.tol) &&
         (Out_data->sweep[0] < Params.sweep_max)){
      i_max = IndMaxAbsDouble(Vars->r, A.n);
      u_max = Vars->u[i_max];
      Vars->u[i_max] += Vars->r[i_max] / A.diag[i_max];
      du = Vars->u[i_max] - u_max;
      for (int i = 0; i < A.j_ptr[i_max+1]-A.j_ptr[i_max]; i++) {
         ind = A.j_ptr[i_max]+i;
         row = A.i[ind];
         elem = A.a[ind];
         Vars->r[row] -= du * elem;
         Out_data->comm[0]++;
         comm++;
      }
      Out_data->relax_hist[i_max]++;
      Out_data->relax[0]++;
      relax++;
      Out_data->sweep[0]++;
   }
   Out_data->relax_scaled[0] += (relax / (double)A.n);
   Out_data->comm_scaled[0] += (comm / (double)A.n);
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
}

/*
 * Sequential Gauss-Seidel
 */
void SEQ_GaussSeidelRes(CSC A, 
                        SolveVars *Vars,
                        SolveParams Params,
                        SolveData *Out_data)
{
   int ind, row;
   double r_norm, u_prev, r_prev, elem, relax, comm;
   double start = omp_get_wtime();
   while((Norm2(Vars->r, A.n) > Params.tol) &&
         (Out_data->sweep[0] < Params.sweep_max)){
      relax = 0; comm = 0;
      for (int i = 0; i < A.n; i++){
         u_prev = Vars->u[i];
         Vars->u[i] += Vars->r[i] / A.diag[i];
         Vars->du[i] = Vars->u[i] - u_prev;
         for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
            ind = A.j_ptr[i]+j;
            row = A.i[ind];
            elem = A.a[ind];
            Vars->r[row] -= Vars->du[i] * elem;
            Out_data->comm[0]++;
            comm++;
         }
         Out_data->relax[0]++;
         relax++;
      }
      Out_data->sweep[0]++;
      Out_data->relax_scaled[0] += (relax / (double)A.n);
      Out_data->comm_scaled[0] += (comm / (double)A.n);
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
}


void SEQ_JacobiRes(CSC A, 
                   SolveVars *Vars,
                   SolveParams Params,
                   SolveData *Out_data)
{
   int ind, row;
   double du_norm, u_prev, elem, relax, du;
   double *r_prev = (double *)malloc(A.n * sizeof(double));
   if (Params.sweep_max > 0){
      Out_data->work_iter[0] = 1;
   }
   else {
      Out_data->work_iter[0] = 0;
   }
   double start = omp_get_wtime();   
   while((Norm2(Vars->r, A.n) > Params.tol) && 
         (Out_data->sweep[0] < Params.sweep_max)){
      relax = 0;
      for (int i = 0; i < A.n; i++) r_prev[i] = Vars->r[i];
      for (int i = 0; i < A.n; i++){
         Vars->u[i] += r_prev[i] / A.diag[i];
         for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
            ind = A.j_ptr[i]+j;
            row = A.i[ind];
            elem = A.a[ind];
            Vars->r[row] -=  (elem * r_prev[i] / A.diag[i]);
         }
         Out_data->relax[0]++;
         relax++;
      }
      Out_data->sweep[0]++;
      Out_data->relax_scaled[0] += (relax / (double)A.n);
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   free(r_prev);
}

void BlockSweepGS(CSC A,
                  OrderInfo P,
                  SolveVars *Vars,
                  SolveData *Out_data,
                  double **r_loc,
                  int p)
{
   int ind, row, shift_row, k;
   double u_prev, elem, start, ra;
   double relax = 0;

   int disp = P.disp[p];
   int part = P.part[p];
   for (int i = P.disp[p]; i < P.disp[p+1]; i++){
      u_prev = Vars->u[i];
      Vars->u[i] += (r_loc[p][i-disp] / A.diag[i]);
      Vars->du[i] = Vars->u[i] - u_prev;
      for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
         ind = A.j_ptr[i]+j;
         row = A.i[ind];
         elem = A.a[ind];
         k = row - disp;
         ra = Vars->du[i] * elem;
         if ((k < part) && (k >= 0)){
            r_loc[p][k] -= ra;
         }
         Vars->r[row] -= ra;
      }
      relax++;
      Out_data->relax[0]++;
   }
   Out_data->relax_scaled[0] += (relax / (double)A.n);
}

void BlockSweepDirect(CSC D,
                      CSC B,
                      PardisoInfo Pard,
                      OrderInfo P,
                      SolveVars *Vars,
                      SolveData *Out_data,
                      double **r_loc,
                      int p)
{
   int ind, row, k;
   double u_prev, elem, start, ra;
   double relax = 0;

   int disp = P.disp[p];
   int part = P.part[p];

   PARDISO(Pard.pt,
           &(Pard.maxfct),
           &(Pard.mnum),
           &(Pard.mtype),
           &(Pard.phase),
           &(Pard.csr.n),
           D.a,
           Pard.csr.ia,
           Pard.csr.ja,
           &(Pard.idum),
          // Pard.perm,
           &(Pard.nrhs),
           Pard.iparm,
           &(Pard.msglvl),
           r_loc[p],
           &(Vars->du[disp]),
           &(Pard.error));
   for (int i = 0; i < B.n; i++){
      u_prev = Vars->u[disp+i];
      Vars->u[disp+i] += Vars->du[disp+i];
      Vars->r[disp+i] = 0;
      relax++;
      Out_data->relax[0]++;
     // for (int j = 0; j < B.j_ptr[i+1]-B.j_ptr[i]; j++) {
     //    ind = B.j_ptr[i]+j;
     //    row = B.i[ind];
     //    elem = B.a[ind];
     //    Vars->r[row] -= (Vars->du[disp+i] * elem);
     // }
   }
   Out_data->relax_scaled[0] += (relax / (double)D.n);
}


void BlockOffDiagDirect(CSC B,
                        OrderInfo P,
                        SolveVars *Vars,
                        int p)
{
   int ind, row, k;
   double u_prev, elem, start, ra;
   double relax = 0;

   int disp = P.disp[p];
   int part = P.part[p];
   for (int i = 0; i < B.n; i++){
      for (int j = 0; j < B.j_ptr[i+1]-B.j_ptr[i]; j++) {
         ind = B.j_ptr[i]+j;
         row = B.i[ind];
         elem = B.a[ind];
         Vars->r[row] -= (Vars->du[disp+i] * elem);
      }
   }
}

void SEQ_BlockJacobiRes(CSC A,
                        CSC *D,
                        CSC *B,
                        PardisoInfo *Pard,
                        OrderInfo P,
                        SolveVars *Vars,
                        SolveParams Params,
                        SolveData *Out_data)
{
   int ind, row;
   double du_norm, u_prev, elem, relax, du;
   double *r_prev = (double *)malloc(A.n * sizeof(double));
   double **r_loc = (double **)calloc(P.nparts, sizeof(double *));

   
   if (Params.sweep_max > 0){
      Out_data->work_iter[0] = 1;
      for (int p = 0; p < P.nparts; p++){
         r_loc[p] = (double *)calloc(P.part[p], sizeof(double)); 
      }
   }
   else {
      free(r_loc);
      Out_data->work_iter[0] = 0;
   }
   
   double start = omp_get_wtime();
   while((Norm2(Vars->r, A.n) > Params.tol) &&
         (Out_data->sweep[0] < Params.sweep_max)){
      relax = 0;
      for (int p = 0; p < P.nparts; p++){
         for (int i = 0; i < P.part[p]; i++){
            r_loc[p][i] = Vars->r[P.disp[p]+i];
         }
      }
      for (int p = 0; p < P.nparts; p++){
         if (loc_direct_flag){
            BlockSweepDirect(D[p], B[p], Pard[p], P, Vars, Out_data, r_loc, p);
         }
         else {
            BlockSweepGS(A, P, Vars, Out_data, r_loc, p);
         }
      }
      for (int p = 0; p < P.nparts; p++){
         BlockOffDiagDirect(B[p], P, Vars, p);
      }
      Out_data->sweep[0]++;
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   for (int p = 0; p < P.nparts; p++){
      free(r_loc[p]);
   }
   free(r_loc);
   free(r_prev);
}

int SEQ_BlockCheckRelax(SEQ_BlockNeighbInfo Neighb, 
                        double *r_norm,
                        int p)
{
   int relax_flag = 1;
   int row, ind, count = 0;
   double r_norm_neighb;
   double r_norm_p = r_norm[p];
   for (int i = 0; i < Neighb.size[p]; i++){
      r_norm_neighb = r_norm[Neighb.blocks[p][i]];
     // printf("%e, %e\n", r_norm_p, r_norm_neighb);
      if (!AlmostEqual2sComplement_DBL(r_norm_neighb, r_norm_p, 1)){
         if (r_norm_neighb > r_norm_p){
            count++;
            if (count == ps_num_relax){
               relax_flag = 0;
               break;
            }
         }
      }
   }
   return relax_flag;
}

void SEQ_BlockParSouthwellRes(CSC A,
                              CSC *D,
                              CSC *B,
                              PardisoInfo *Pard,
                              OrderInfo P,
                              SolveVars *Vars,
                              SolveParams Params,
                              SolveData *Out_data)
{
   int ind, row, relax_flag;
   double du_norm, u_prev, elem, relax, du;
   double *r_prev = (double *)malloc(A.n * sizeof(double));
   double **r_loc = (double **)calloc(P.nparts, sizeof(double *));
   double *r_norm = (double *)calloc(P.nparts, sizeof(double *));
   int *relax_mask = (int *)calloc(P.nparts, sizeof(int *));
   SEQ_BlockNeighbInfo Neighb;

   SEQ_SetupNeighbBlocks(P, B, &Neighb);

   
   if (Params.sweep_max > 0){
      for (int p = 0; p < P.nparts; p++){
         r_loc[p] = (double *)calloc(P.part[p], sizeof(double)); 
      }
   }
   else {
      free(r_loc);
      Out_data->work_iter[0] = 0;
   }
   
   double start = omp_get_wtime();
   while((Norm2(Vars->r, A.n) > Params.tol) &&
         (Out_data->sweep[0] < Params.sweep_max)){
      relax = 0;
      for (int p = 0; p < P.nparts; p++){
         for (int i = 0; i < P.part[p]; i++){
            r_loc[p][i] = Vars->r[P.disp[p]+i];
            if (loc_direct_flag && relax_mask[p]){
               Vars->r[P.disp[p]+i] = 0;
            }
         }
         r_norm[p] = Norm2(r_loc[p], P.part[p]);
      }
      for (int p = 0; p < P.nparts; p++){
         relax_mask[p] = 0;
         relax_flag = SEQ_BlockCheckRelax(Neighb, r_norm, p);
         if (relax_flag){
            relax_mask[p] = 1;
            if (loc_direct_flag){
               BlockSweepDirect(D[p], B[p], Pard[p], 
                                P, Vars, Out_data, r_loc, p);
            }
            else {
               BlockSweepGS(A, P, Vars, Out_data, r_loc, p);
            }
         }
      }
      for (int p = 0; p < P.nparts; p++){
         if (relax_mask[p]){
            BlockOffDiagDirect(B[p], P, Vars, p);
         }
      }
      Out_data->sweep[0]++;
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   for (int p = 0; p < P.nparts; p++){
      free(r_loc[p]);
   }
   Out_data->work_iter[0] =
      ((double)SumInt(relax_mask, P.nparts))/(double)P.nparts;
   free(r_norm);
   free(r_loc);
   free(r_prev);
}

void SMEM_JacobiRes(CSC A,
                    SolveVars *Vars,
                    SolveParams Params,
                    SolveData *Out_data)
{ 
   if (threads == 1){
      SEQ_JacobiRes(A, Vars, Params, Out_data);
      return;
   }

   double *r_prev = (double *)calloc(A.n, sizeof(double));

   double start = omp_get_wtime();
   #pragma omp parallel num_threads(threads)
   {
      int ind, row, relax_flag;
      double relax, u_prev, elem;
      while((OMP_Norm2(Vars->r, A.n) > Params.tol) && 
            (Out_data->sweep[0] < Params.sweep_max)){
         #pragma omp for 
         for (int i = 0; i < A.n; i++){
            r_prev[i] = Vars->r[i];
         }
         relax = 0;
         #pragma omp for
         for (int i = 0; i < A.n; i++){
            Vars->u[i] += r_prev[i] / A.diag[i];
            for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
               ind = A.j_ptr[i]+j;
               row = A.i[ind];
               elem = A.a[ind];
               #pragma omp atomic
               Vars->r[row] -= (elem * r_prev[i] / A.diag[i]);
            }
            #pragma omp atomic
            Out_data->relax[0]++;
            relax++;
         }
         #pragma omp single
         Out_data->sweep[0]++;
         #pragma omp atomic
         Out_data->relax_scaled[0] += (relax / (double)A.n);
      }
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   free(r_prev);
}

int SEQ_CheckRelax(CSC A, double *r_prev, int i)
{
   int relax_flag = 1;
   int row, ind, q = 1;
   double r_row;
   double r_i = fabs(r_prev[i]);
   for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++){
      ind = A.j_ptr[i]+j;
      row = A.i[ind];
      if (row != i){
         r_row = fabs(r_prev[row]);
         if (AlmostEqual2sComplement_DBL(r_row, r_i, 1)){
            if (i < row){
               relax_flag = 0;
               break;
            }
         }
         else {
            if (r_row > r_i){
               relax_flag = 0;
               break;
            }
         }
      }
   }
   return relax_flag;
}

void SEQ_ParSouthwellRes(CSC A,
                         SolveVars *Vars,
                         SolveParams Params,
                         SolveData *Out_data)
{
   int ind, row, max_degree, min_degree;
   double relax, comm, u_prev, elem;
 
   double *r_prev = (double *)calloc(A.n, sizeof(double));

   int relax_flag;
   double start = omp_get_wtime();
   Out_data->work_iter[0] = 0;
   while((Norm2(Vars->r, A.n) > Params.tol) && 
         (Out_data->sweep[0] < Params.sweep_max)){
      for (int i = 0; i < A.n; i++){
         r_prev[i] = Vars->r[i];
         Out_data->relax_mask[i] = 0;
         Out_data->relax_hist[i] = 0;
      }
      relax = 0;
      for (int i = 0; i < A.n; i++){
         Out_data->relax_mask[i] = 0;
         relax_flag = SEQ_CheckRelax(A, r_prev, i);
         if (relax_flag){
            Vars->u[i] += r_prev[i] / A.diag[i];
            for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
               ind = A.j_ptr[i]+j;
               row = A.i[ind];
               elem = A.a[ind];
               Vars->r[row] -=  (elem * r_prev[i] / A.diag[i]);
            }
            Out_data->relax_mask[i] = 1;
            Out_data->relax_hist[i]++;
            Out_data->relax[0]++;
            relax++;
         } 
      }
      Out_data->sweep[0]++;
      Out_data->relax_scaled[0] += (relax / (double)A.n);
   }
   Out_data->work_iter[0] =
      ((double)SumInt(Out_data->relax_mask, A.n))/(double)A.n;
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   free(r_prev);
}


void SMEM_ParSouthwellRes(CSC A,
                          SolveVars *Vars,
                          SolveParams Params,
                          SolveData *Out_data)
{ 
   if (threads == 1){
      SEQ_ParSouthwellRes(A, Vars, Params, Out_data);
      return;
   }

   double *r_prev = (double *)calloc(A.n, sizeof(double));

   double start = omp_get_wtime();
   #pragma omp parallel num_threads(threads)
   {
      int ind, row, relax_flag;
      double u_prev, elem, relax;
      while((OMP_Norm2(Vars->r, A.n) > Params.tol) && 
            (Out_data->sweep[0] < Params.sweep_max)){
         #pragma omp for 
         for (int i = 0; i < A.n; i++){
            r_prev[i] = Vars->r[i];
            Out_data->relax_mask[i] = 0;
         }
         relax = 0;
         #pragma omp for
         for (int i = 0; i < A.n; i++){
            relax_flag = SEQ_CheckRelax(A, r_prev, i);
            if (relax_flag){
               Vars->u[i] += r_prev[i] / A.diag[i];
               for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
                  ind = A.j_ptr[i]+j;
                  row = A.i[ind];
                  elem = A.a[ind];
                  #pragma omp atomic
                  Vars->r[row] -= (elem * r_prev[i] / A.diag[i]);
               }
               Out_data->relax_mask[i] = 1;
               Out_data->relax_hist[i]++;
               #pragma omp atomic
               Out_data->relax[0]++;
               relax++;
            } 
         }
         #pragma omp single
         Out_data->sweep[0]++;
         #pragma omp atomic
         Out_data->relax_scaled[0] += (relax / (double)A.n);
      }
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
   free(r_prev);
}

void SEQ_MCGS(CSC A,
              OrderInfo P,
              SolveVars *Vars,
              SolveParams Params,
              SolveData *Out_data)

{
   int row, disp, ind;
   int break_flag = 0;
   double relax;
   double u_prev, elem;
   Out_data->work_iter[0] = 0;
   double start = omp_get_wtime();
   while(1){
      for (int c = 0; c < P.nparts; c++){
         if ((Norm2(Vars->r, A.n) < Params.tol) ||
             (Out_data->sweep[0] >= Params.sweep_max)){
            break_flag = 1;
            break;
         }
         relax = 0;
         for (int i = P.disp[c]; i < P.disp[c]+P.part[c]; i++){
            Out_data->work_iter[0] = (double)P.part[c]/(double)A.n;
            u_prev = Vars->u[i];
            Vars->u[i] += Vars->r[i] / A.diag[i];
            Vars->du[i] = Vars->u[i] - u_prev;
            for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
               ind = A.j_ptr[i]+j;
               row = A.i[ind];
               elem = A.a[ind];
               Vars->r[row] -= Vars->du[i] * elem;
            }
            Out_data->relax[0]++;
            relax++;
         }
         Out_data->sweep[0]++;
         Out_data->relax_scaled[0] += (relax / (double)A.n);
      }
      if (break_flag) break;
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
}

void SMEM_MCGS(CSC A,
               OrderInfo P,
               SolveVars *Vars,
               SolveParams Params,
               SolveData *Out_data)

{
   if (threads == 1){
      SEQ_MCGS(A, P, Vars, Params, Out_data); 
      return;
   }
   double start = omp_get_wtime();
   #pragma omp parallel num_threads(threads)
   {
      int row, disp, ind;
      double relax;
      double u_prev, elem;
      while((OMP_Norm2(Vars->r, A.n) > Params.tol) &&
            (Out_data->sweep[0] < Params.sweep_max)){
         for (int c = 0; c < P.nparts; c++){
            relax = 0;
            #pragma omp for
            for (int i = P.disp[c]; i < P.disp[c]+P.part[c]; i++){
               u_prev = Vars->u[i];
               Vars->u[i] += Vars->r[i] / A.diag[i];
               Vars->du[i] = Vars->u[i] - u_prev;
               for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++) {
                  ind = A.j_ptr[i]+j;
                  row = A.i[ind];
                  elem = A.a[ind];
                  #pragma omp atomic
                  Vars->r[row] -= Vars->du[i] * elem;
               }
               #pragma omp atomic
               Out_data->relax[0]++;
               relax++;
            }
            #pragma omp single
            Out_data->sweep[0]++;
            #pragma omp atomic
            Out_data->relax_scaled[0] += (relax / (double)A.n);
         } 
      }
   }
   Out_data->wtime_tot[0] = omp_get_wtime() - start;
}
