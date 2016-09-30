#include <stdio.h>
#include <stdlib.h>

#define NTHRDS 4
#include <omp.h>

static long num_p = 100000; // Numero total de segmentos
double paso;               // Delta x

int main(void){
   int i;
   double x, pi, suma = 0.0;

   paso = 1.0 / (double)num_p;

   // Se calcula la suma
   omp_set_num_threads(NTHRDS);
#pragma omp parallel for private(x) reduction(+:suma)
   for(i = 0; i < num_p; i++){
      x = (i + 0.5) * paso;
      suma = suma + 4.0 / (1.0 + x * x);
   }

   // Se estima pi
   pi = paso * suma;

  printf("El valor aproximado de pi es: %3.10lf \n", pi);

   return EXIT_SUCCESS;
}

