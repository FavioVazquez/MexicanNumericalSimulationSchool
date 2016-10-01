#include <stdio.h>
#include <stdlib.h>

#define NTHRDS 4
#include <omp.h>

static long num_p = 10000; // Numero total de segmentos
double paso;               // Delta x

int main(void){
   int num_t;
   double pi, suma[4], sumaf;

   paso = 1.0 / (double)num_p;

   // Se calcula la suma en paralelo
#pragma omp parallel num_threads(NTHRDS)
{
   int i;
   int mi_id = omp_get_thread_num();
   double x;

   num_t = omp_get_num_threads();
   suma[mi_id]=0.0;

   for(i = 0 + mi_id; i < num_p; i=i+num_t){
      x = (i + 0.5) * paso;
      suma[mi_id] = suma[mi_id] + 4.0 / (1.0 + x * x);
   }
}

   // Se estima pi de forma secuencial
   int i;
   for(i=0;i<NTHRDS;i++)
      sumaf=sumaf+suma[i];

   pi = paso * sumaf;

  printf("El valor aproximado de pi es: %3.10lf \n", pi);

   return EXIT_SUCCESS;
}

