#include "Southwell.h"

void List_to_Metis(MetisGraph *G,
                   Triplet *T,
                   std::vector<std::list<int>> col_list,
                   std::vector<std::list<double>> elem_list);

//void ReadText_fscanf_metis(FILE *mat_file, MetisGraph *G)
//{
//   int ip = 0, ip_i, i, j, row, col, insert_flag, file_len = 0;
//   int col_i, k, n;
//   double elem;
//   List *col_list;
//   int max_row = 0; 
//   rewind(mat_file);
//   G->nnz = 0;
//   while(fscanf(mat_file, "%d %d %lg", &row, &col, &elem) == 3){
//      if (elem != 0.){
//         file_len++;
//         G->nnz++;
//         if (row > max_row) max_row = row;
//      }
//   }
//   if (file_len <= 0){
//      printf("ReadMatrix() Error: file is empty.\n");
//      return;
//   }
//   G->n = (idx_t)max_row;
//   n = row;
//
//   col_list = (List *)malloc(n * sizeof(List));
//   for (int i = 0; i < n; i++){ 
//      col_list[i] = newList();
//      moveTo(col_list[i], 0);
//   }
//
//   G->xadj = (idx_t *)malloc((n+1) * sizeof(idx_t)); 
//   rewind(mat_file);
//   while(fscanf(mat_file, "%d %d %lg", &row, &col, &elem) == 3){
//      if (elem != 0.){
//         append(col_list[col-1], row-1);
//         moveNext(col_list[col-1]);
//      }
//   }
//   G->adjncy = (idx_t *)malloc((int)G->nnz * sizeof(idx_t));
//   G->adjwgt = (real_t *)malloc((int)G->nnz * sizeof(real_t));
//   k = 0;
//   G->xadj[0] = 0;
//   for (int i = 0; i < n; i++){
//       moveTo(col_list[i], 0);
//       G->xadj[i+1] = G->xadj[i] + (idx_t)length(col_list[i]);
//       for (int j = 0; j < length(col_list[i]); j++){
//          G->adjncy[k] = (idx_t)getElement(col_list[i]);
//          moveNext(col_list[i]);
//          k++;
//       }
//   }
//   G->xadj[n] = G->nnz;
//
//   int *row_count = (int *)calloc(n, sizeof(int));
//   rewind(mat_file);
//   while(fscanf(mat_file, "%d %d %lg", &row, &col, &elem) == 3){
//      if (elem != 0.){
//         k = row_count[col-1];
//         j = G->xadj[col-1];
//         G->adjwgt[j+k] = elem;
//         row_count[col-1]++;
//      }
//   }
//
//   for (int i = 0; i < n; i++){
//      freeList(&col_list[i]);
//   }
//   free(col_list);
//   free(row_count);
//}

void get_tokens(size_t count, char *str, char *delim, char ***out)
{
   size_t idx  = 0;
   char *token = strtok(str, delim);

   while (token)
   {
       assert(idx < count);
       *(*out + idx++) = strdup(token);
       token = strtok(0, delim);
   }
   assert(idx == count - 1);
   *(*out + idx) = 0;
}

size_t string_split(char *str, const char a_delim, char ***out)
{
   size_t count     = 0;
   char *tmp        = str;
   char *prev_delim = NULL;
   char delim[2];
   delim[0] = a_delim;
   delim[1] = 0;

   /* Count how many elements will be extracted. */
   while (*tmp)
   {
       if (a_delim == *tmp)
       {
           count++;
           prev_delim = tmp;
       }
       tmp++;
   }

   /* Add space for trailing token. */
   count += prev_delim < (str + strlen(str) - 1);

   /* Add space for terminating null string so caller
    * knows where the list of returned strings ends. */
   count++;

   *out = (char **)malloc(count * sizeof(char *));

   if (*out)
   {
      get_tokens(count, str, delim, out);
   }

   return count-1;
}

/* INCOMPLETE FUNCTION */
void ReadText_fread_metis(FILE *mat_file, MetisGraph *G, Triplet *T)
{
   using namespace std;
   int col_max = 0;
   int temp;
   size_t result, size;
   size_t num_lines, num_elems;
   char *buffer_char;
   char **elem_tokens, **line_tokens;

   fseek(mat_file , 0 , SEEK_END);
   size = ftell(mat_file);
   rewind(mat_file);
   buffer_char = (char *)malloc(sizeof(char) * size);
   fread(buffer_char, sizeof(char), size, mat_file);
   num_lines = string_split(buffer_char, '\n', &line_tokens);
   T->i = (int *)malloc(num_lines * sizeof(int));
   T->j = (int *)malloc(num_lines * sizeof(int));
   T->a = (double *)malloc(num_lines * sizeof(double));

   T->nnz = 0;
   for (int k = 0; k < num_lines; k++){
      num_elems = string_split(line_tokens[k], ' ', &elem_tokens);
      if (num_elems == 3){
         T->i[k] = atoi(elem_tokens[0]);
         T->j[k] = atoi(elem_tokens[1]);
         T->a[k] = atof(elem_tokens[2]);
         T->nnz++;
      }
      else break;
      for (int i = 0; i < num_elems; i++) free(elem_tokens[i]);
      free(elem_tokens);
   }

   T->n = *max_element(T->i, T->i+T->nnz); 

   printf("%d, %d\n", T->nnz, T->n);
   
   for (int i = 0; i < num_lines; i++) free(line_tokens[i]);
   free(line_tokens);
}

void ReadBinary_fread_metis(FILE *mat_file, MetisGraph *G, Triplet *T)
{
   using namespace std;
   size_t size;
   int temp_size;
   int k, q;
   int row, col;
   double elem;
   Triplet_AOS *buffer;

   fseek(mat_file , 0 , SEEK_END);
   size = ftell(mat_file);
   rewind(mat_file);
   buffer = (Triplet_AOS *)malloc(sizeof(Triplet_AOS) * size);
   fread(buffer, sizeof(Triplet_AOS), size, mat_file);

   
   T->nnz = G->nnz = size/sizeof(Triplet_AOS);
   int *row_array = (int *)malloc(T->nnz * sizeof(int));
   std::vector<list<int>> col_list(G->nnz);
   std::vector<list<double>> elem_list(G->nnz);
   for (int k = 0; k < T->nnz; k++){
      row = buffer[k].i;
      col = buffer[k].j;
      elem = buffer[k].a;
      row_array[k] = row;
      col_list[col-1].push_back(row-1);
      elem_list[col-1].push_back(elem);
   }
   G->n = T->n = *max_element(row_array, row_array+T->nnz);

   G->xadj = (idx_t *)malloc(((int)G->n+1) * sizeof(idx_t));
   G->adjncy = (idx_t *)malloc((int)G->nnz * sizeof(idx_t));
   G->adjwgt = (real_t *)malloc((int)G->nnz * sizeof(real_t));
  
   T->i = (int *)malloc(T->nnz * sizeof(int));
   T->j = (int *)malloc(T->nnz * sizeof(int));
   T->a = (double *)malloc(T->nnz * sizeof(double));

   List_to_Metis(G, T, col_list, elem_list);  
}

void Write_csc(CSC A, int base)
{
   int row, col, k;
   double elem;
   char buffer[100];
   if (color_flag){
      strcpy(buffer, "mc_matrix_matlab.txt");
   }
   else {
      strcpy(buffer, "metis_matrix_matlab.txt");
   }
   FILE *out_file;
   remove(buffer);
   out_file = fopen(buffer, "a");
   for (int i = 0; i < A.n; i++){
      k = A.j_ptr[i];
      for (int j = 0; j < A.j_ptr[i+1]-A.j_ptr[i]; j++){
         row = A.i[k+j];
         col = i;
         elem = A.a[k+j];
         fprintf(out_file, "%d   %d   %e\n", col+base, row+base, elem);
      }
   }
   fclose(out_file);
}

void WriteBlocks_csc(char *buffer,
                     CSC *D,
                     CSC *B,
                     OrderInfo P,
                     int base)
{
   int row, col, k;
   double elem;
   FILE *out_file;
   out_file = fopen(buffer, "w");
   for (int p = 0; p < P.nparts; p++){
      for (int i = 0; i < D[p].n; i++){
         k = D[p].j_ptr[i];
         for (int j = 0; j < D[p].j_ptr[i+1]-D[p].j_ptr[i]; j++){
            row = D[p].i[k+j];
            col = i;
            elem = D[p].a[k+j];
            fprintf(out_file, "%d %d %e\n",
                    row+base, P.disp[p]+col+base, elem);
         }
      }
      for (int i = 0; i < B[p].n; i++){
         k = B[p].j_ptr[i];
         for (int j = 0; j < B[p].j_ptr[i+1]-B[p].j_ptr[i]; j++){
            row = B[p].i[k+j];
            col = i;
            elem = B[p].a[k+j];
            fprintf(out_file, "%d %d %e\n",
                    row+base, P.disp[p]+col+base, elem);
         }
      }
   }
   fclose(out_file);
}

void Write_triplet(Triplet T, int base)
{
   int row, col, k, rank_p = 0, num_p = 1;
   double elem;
   char buffer[100];
   if (color_flag){
      strcpy(buffer, "mc_matrix_matlab.txt");
   }
   else {
      strcpy(buffer, "metis_matrix_matlab.txt");
   }
   FILE *out_file;
   remove(buffer);
   out_file = fopen(buffer, "a");
   for (int i = 0; i < T.nnz; i++){
      fprintf(out_file, "%d   %d   %e\n", T.i[i]+base, T.j[i]+base, T.a[i]);
   }
   fclose(out_file);
}

void Write_metis(FILE *out_file, MetisGraph G, int base)
{
   int row, col, k, rank_p = 0;
   double elem;
   int q = 0;
   for (int i = 0; i < G.n; i++){
      k = G.xadj[i];
      for (int j = 0; j < G.xadj[i+1]-G.xadj[i]; j++){
         row = G.adjncy[k+j];
         col = i;
         elem = G.adjwgt[k+j];
         fprintf(out_file, "%d   %d   %e\n", col+base, row+base, elem);
         q++;
      }
   }
}

void Print_csc(CSC A)
{
   for (int i = 0; i < A.nnz; i++){
      printf("%d, %.5f\n", i, A.a[i]);
   }
}

void Print_metis(MetisGraph G)
{
   for (int i = 0; i < (int)G.nnz; i++){
      printf("%d, %d\n", i, G.adjwgt[i]);
   }
}

void ReadText_fscanf_array(char *file_name, double **v, int n)
{
   double elem;
   int n_test = 0;
   FILE *file = fopen(file_name, "r");
   while(fscanf(file, "%lg", &elem) == 1){
      n_test++;
   }
   if (n_test != n){
      printf("Vector read error: std::vector size passed to ReadVector() is different than actual size\n"
             "size is %d, but should be %d\n", n_test, n);
      exit(1);
   }
   rewind(file);
   int i = 0;
   while(fscanf(file, "%lg", &elem) == 1){
      (*v)[i] = elem;
      i++;
   }
}

void List_to_Metis(MetisGraph *G,
                   Triplet *T,
                   std::vector<std::list<int>> col_list,
                   std::vector<std::list<double>> elem_list)
{
   int temp_size;
   int k = 0;
   G->xadj[0] = 0;
   for (int i = 0; i < G->n; i++){
      col_list[i].begin();
      elem_list[i].begin();
      temp_size = col_list[i].size();
      G->xadj[i+1] = G->xadj[i] + (idx_t)temp_size;
      for (int j = 0; j < temp_size; j++){
         G->adjncy[k] = col_list[i].front();
         col_list[i].pop_front();
         G->adjwgt[k] = elem_list[i].front();
         elem_list[i].pop_front();

         T->j[k] = i;
         T->i[k] = G->adjncy[k];
         T->a[k] = G->adjwgt[k];

         k++;
      }
   }
}

void ReadVector(char *file_name, double **v, int n)
{
   double elem;
   int n_test = 0;
   FILE *file = fopen(file_name, "r");
   while(fscanf(file, "%lg", &elem) == 1){
      n_test++;
   }
   if (n_test != n){
      printf("Vector read error: vector size passed to ReadVector() is different than actual size\n"
             "size is %d, but should be %d\n", n_test, n);
      exit(1);
   }
   rewind(file);
   int i = 0;
   while(fscanf(file, "%lg", &elem) == 1){
      (*v)[i] = elem;
      i++;
   }
}
