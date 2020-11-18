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
void gen_data_arr(double * array, int size, int real_size);
void gen_data_vec(double * array, int size, int real_size);
/* función para multiplicar iterativamente un matriz 
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double* local_A, double* x, double* y, double* local_y, int n, int it, int rows_per_proc, int rank);
/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char* name, double*  y, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;
  int n, real_n, iters,rows_per_proc;
  long seed;

  double* local_A=NULL;
  double* local_y=NULL;

  int my_rank,comm_sz;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  // Obtener las dimensiones
  if(my_rank==0){
    printf("Ingrese la dimensión n:\n");
    scanf("%d", &n);
    printf("Ingrese el número de iteraciones:\n");
    scanf("%d", &iters);
    printf("Ingrese semilla para el generador de números aleatorios:\n");
    scanf("%ld", &seed);
    real_n=n;
  }

  MPI_Bcast ( &n, 1 , MPI_INT , 0 , MPI_COMM_WORLD );
  MPI_Bcast ( &iters, 1 , MPI_INT , 0 , MPI_COMM_WORLD );
  
  // se hace n divisible por el numero de procesos
  if(n%comm_sz!=0){
      n=n+(comm_sz-(n%comm_sz));
  }
  x = malloc(sizeof(double) * n);
  y = malloc(sizeof(double) * n);
  if(my_rank==0){
    srand(seed); 
    A = malloc(sizeof(double) * n * n);
    gen_data_arr(A, n, real_n);
    gen_data_vec(x, n, real_n);
  }

  MPI_Bcast(x, n , MPI_INT , 0 , MPI_COMM_WORLD );
  rows_per_proc=n/comm_sz;
  local_A=malloc(sizeof(double)*rows_per_proc*n);
  local_y=malloc(sizeof(double)*rows_per_proc);
  MPI_Scatter(A, rows_per_proc*n, MPI_DOUBLE, local_A, rows_per_proc*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

  //generar valores para las matrices
 
  if(my_rank==0){
    //print_vector("x", x, n);
    //print_vector("A", A, n*n);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  mat_vect_mult(local_A, x, y, local_y, n, iters, rows_per_proc, my_rank);

  MPI_Gather( local_y , rows_per_proc, MPI_DOUBLE , y , rows_per_proc, MPI_DOUBLE , 0, MPI_COMM_WORLD );
  if(my_rank==0){
    print_vector("y", y, real_n);
    free(A);
  }
  
  free(x);
  free(y);
  MPI_Finalize();
  return 0;
}

void gen_data_arr(double * array, int size, int real_size){
  int i,j;
  for (i = 0; i < size; i++){
      for (j=0; j<size; j++){
          if(i<real_size && j<real_size){
              array[i*size+j] = (double) rand() / (double) RAND_MAX;
          }
          else{
              array[i*size+j] = 0;
          }
      }
  }
}

void gen_data_vec(double * array, int size, int real_size){
  int i;
  for (i = 0; i < size; i++){
    if(i<real_size){
        array[i] = (double) rand() / (double) RAND_MAX;
    }
    else{
        array[i] = 0;
    }
  }
}

void mat_vect_mult(double* local_A, double* x, double* y, double* local_y, int n, int it, int rows_per_proc, int rank){
  int h, i, j;
  for(h = 0; h < it; h++){ 
    for(i=0; i<rows_per_proc; i++){
        local_y[i]=0.0;
        for(j=0; j<n; j++){
            local_y[i]+=local_A[i*n+j]*x[j];
        }
    }
    if(rank==1){
        print_vector("local_A",local_A, rows_per_proc*n);
        print_vector("x",x,n);
        print_vector("local_y",local_y,rows_per_proc);
    }
    MPI_Allgather( local_y , rows_per_proc , MPI_DOUBLE , y , rows_per_proc , MPI_DOUBLE , MPI_COMM_WORLD );
    // x <= y
    for(i = 0; i < n; i++)
      x[i] = y[i];
  }
}

void print_vector(char* name, double*  y, int m) {
   int i;
   printf("\nVector %s\n", name);
   for (i = 0; i < m; i++)
      printf("%f ", y[i]);
   printf("\n");
}
