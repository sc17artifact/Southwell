#include "Southwell.h"
#include "Misc.h"
#include "pseudo_peripheral_node.h"

extern int threads;
extern int mat_file_flag;
extern int color_flag;
extern int reorder_flag;
extern int format_out_flag;

int BreadthFirstColor(CSC A, OrderInfo *P, int start);
int BalanceColorParts(CSC A, OrderInfo *P);

void Multicolor(CSC A, OrderInfo *P)
{
   int start = 1;
   int level_num;
  // int *level = (int *)calloc(A.n, sizeof(int));
  // int *level_row = (int *)calloc(A.n+1, sizeof(int));
  // int *mask = (int *)calloc(A.n, sizeof(int));
  // for (int i = 0; i < A.n; i++) mask[i] = 1;
  // root_find(&start, A.nnz, A.j_ptr, A.i, 
  //           mask, &level_num, level_row, level, A.n);
  // start--;
   start = 0;
   BreadthFirstColor(A, P, start);
   P->part = (int *)calloc(P->nparts, sizeof(int));
   for (int i = 0; i < A.n; i++){
      P->part[P->perm[i]]++;
   }
   BalanceColorParts(A, P);

  // free(level);
  // free(level_row);
  // free(mask);
}

int BalanceColorParts(CSC A, OrderInfo *P)
{  
   int min_color, found_color, row; 
   int *colors = (int *)malloc(P->nparts * sizeof(int));
   int *temp_part = (int *)malloc(P->nparts * sizeof(int));
   for (int i = 0; i < P->nparts; i++) colors[i] = i;
   for (int i = 0; i < A.n; i++){
     // printf("\nBEFORE\n");
     // for (int j = 0; j < P->nparts; j++) printf("%d, %d\n", P->part[j], colors[j]);
      for (int j = 0; j < P->nparts; j++){
         colors[j] = j;
         temp_part[j] = P->part[j];
      }
      QuicksortPair_int_int(temp_part, colors, 0, P->nparts-1);
     // printf("\nAFTER\n");
     // for (int j = 0; j < P->nparts; j++) printf("%d, %d\n", P->part[j], colors[j]);
     // if (i == 100) break;
      if (P->perm[i] != colors[0]){
         P->part[P->perm[i]]--;
         for (int j = 0; j < P->nparts; j++){
            P->perm[i] = colors[j];
            found_color = 1;
            for (int k = A.j_ptr[i]; k < A.j_ptr[i+1]; k++){
               row = A.i[k];
               if ((P->perm[i] == P->perm[row]) && (i != row)){
                  found_color = 0;
                  break; 
               }
            }
            if (found_color) break;
         }
        // if (!found_color) printf("%d, %d\n", P->perm[i], found_color);
         P->part[P->perm[i]]++;
      }
   }
  // for (int i = 0; i < P->nparts; i++) printf("%d\n", P->part[i]);
   free(colors);
   free(temp_part);
}

int BreadthFirstColor(CSC A, OrderInfo *P, int start)
{
    int i, j, beg, end, nextbeg, nextend, ncol;
    int *visited;
    int nlev, curr_col;
    int n = A.n;
    int *levptr = (int *) malloc((n+1)*sizeof(int));

    P->nparts = 1;
    curr_col = 0;
    nlev = 0;
    levptr[0] = 0;

    // queue
    int *p = (int *) malloc(n*sizeof(int));

    // initialize colors
    for (i=0; i<n; i++)
        P->perm[i] = -1;

    // process first element
    int first = start; // arbitrary
    beg = 0;
    end = 0;
    p[end] = first;
    P->perm[first] = 0;
    end++;

    while (beg < end)
    {
        nlev++;
        levptr[nlev] = end;

        nextbeg = end;
        nextend = end;

        // loop over elements in current level (already in p)
        for (i=beg; i<end; i++)
        {
            int pind = p[i];
            for (j=A.j_ptr[pind]; j<A.j_ptr[pind+1]; j++)
            {
                int ind = A.i[j];

                // ignore any diagonal entry
                if (pind == ind)
                    continue;

                // if not already visited
                if (P->perm[ind] == -1)
                {
                    p[nextend] = ind;
                    nextend++;

                    // propose color
                    // kludge on max number of colors
                    for (int col=0; col<10000; col++)
                    {
                        int jj;
                        int good = 1;
                        for (jj=A.j_ptr[ind]; jj<A.j_ptr[ind+1]; jj++)
                        {
                            if (P->perm[A.i[jj]] == (idx_t)col)
                            {
                                good = 0;
                                break;
                            }
                        }
                        if (good)
                        {
                            if (col > curr_col){
                               P->nparts++;
                               curr_col = col;
                            }
                            P->perm[ind] = (idx_t)col;
                            break;
                        }
                    }
                }
            }
        }
        beg = nextbeg;
        end = nextend;
    }

    free(p);
    free(levptr);

    // check that everything was processed (graph was connected)
    if (end != n)
    {
        printf("number of levels %d\n", nlev);
        printf("graph not connected: end: %d for n %d\n", end, n);
        return -1;
    }

    return 0;
}
