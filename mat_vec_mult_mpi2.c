/**
 *   \file my_it_mat_vect_mult.c
 *   \brief Multiplica iterativamente un matriz nxn 
 *          por un vector de n posiciones
 *
 *   \author Danny Múnera
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double * array, int size);
/* función para multiplicar iterativamente un matriz 
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double* A, double* x, double* y, int n, int it, int p);
/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char* name, double*  y, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;

  double* local_x=NULL;
  double* local_y=NULL;

  int n,local_n, iters;
  int pid, p;
  long seed;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if(pid == 0){
  // Obtener las dimensiones
  printf("Ingrese la dimensión n:\n");
  scanf("%d", &n);
  printf("Ingrese el número de iteraciones:\n");
  scanf("%d", &iters);
  printf("Ingrese semilla para el generador de números aleatorios:\n");
  scanf("%ld", &seed);
  x = malloc(sizeof(double) * n);
  y = malloc(sizeof(double) * n);
  srand(seed);
  gen_data(x, n);
  }

  MPI_Bcast(&seed,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast ( &n, 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
  srand(seed);
  // la matriz A tendrá una representación unidimensional
  A = malloc(sizeof(double) * n * n);
  gen_data(A, n*n);

  
  
  //generar valores para las matrices
  
  if(pid==0){
    //print_vector("x", x, n);
    //print_vector("A", A, n*n);
  }

 

  

  // Local data
  local_n=n/p;
  local_x = malloc(sizeof(double) * local_n);
  local_y = malloc(sizeof(double) * local_n);

  

  MPI_Scatter(x , local_n , MPI_DOUBLE , local_x , local_n , MPI_DOUBLE , 0, MPI_COMM_WORLD );

 

  mat_vect_mult(A, local_x, local_y, n, iters,p);

  MPI_Barrier(MPI_COMM_WORLD);
  if(pid==0){
      printf("llega despues de la multi");
  }

  MPI_Finalize();
  return 0;

  MPI_Gather( local_y , n/p , MPI_DOUBLE , y , n/p , MPI_DOUBLE , 0, MPI_COMM_WORLD );

  MPI_Finalize();
  print_vector("y", y, n);
  free(A);
  if(pid==0){
    free(x);
    free(y);
  }
  
  return 0;
}

void gen_data(double * array, int size){
  int i;
  for (i = 0; i < size; i++)
    array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double* A, double* local_x, double* local_y, int n, int it, int p){
  printf("----entra----\n");
  int h, i, j,k,local_n;
  local_n=n/p;
  for(h = 0; h < it; h++){

      for(i = 0; i < local_n; i++){
        local_y[i] = 0.0;
        for(j = 0; j < local_n; j++){
            local_y[i] += A[i*n+j] * local_x[j];
        }
      }
	      
    /*
    for (i=0 ; i< (n/p) ; i++ )
    {
        for (j = 0 ; j<n ; j++ )
        {
            local_y[i*n+j] = 0.0 ;
            for ( int k = 0 ; k<n ; k++ )
            {
                local_y[i*n+j] += local_x[i*n+k]*A[k*n+j];
            }
        }
    }
    */
    // x <= y
    for(i = 0; i < local_n; i++){
        local_x[i] = local_y[i];
    }
  }
  printf("----sale----\n");
}

void print_vector(char* name, double*  y, int m) {
   int i;
   printf("\nVector %s\n", name);
   for (i = 0; i < m; i++)
      printf("%f ", y[i]);
   printf("\n");
}
