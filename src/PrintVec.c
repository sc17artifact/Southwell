#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <mkl.h>

void RandDouble(double *v,
                int n,
                double low,
                double high)
{
   VSLStreamStatePtr stream;
   vslNewStream(&stream, VSL_BRNG_SFMT19937, 0);
   vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, v, low, high);
}

int main(int argc, char *argv[])
{
   double low = -1, high = 1, rand_range = high-low;
   int ones_flag = 0;
   int rand_flag = 0;
   int gauss_flag = 0;
   int n = 10;
   int arg_iter = 0;
   char file_str[100] = "vector.txt";

   while (arg_iter < argc){
      if (strcmp(argv[arg_iter], "-n") == 0){
         arg_iter++;
         n = atoi(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-rand_range") == 0){
         arg_iter++;
         low = atof(argv[arg_iter++]);
         high = atof(argv[arg_iter]);
      }
      else if (strcmp(argv[arg_iter], "-out") == 0){
         arg_iter++;
         strcpy(file_str, argv[arg_iter]);;
      }
      arg_iter++;
   }

   FILE *file = fopen(file_str, "w");
   double *v = (double *)calloc(n, sizeof(double));
   RandDouble(v, n, low, high);
   
   for (int i = 0; i < n; i++) fprintf(file, "%e\n", v[i]);
    
   fclose(file); 
   return 0;
}
