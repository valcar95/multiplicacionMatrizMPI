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
#include <string.h>

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double * array, int size);
/* función para multiplicar iterativamente un matriz 
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double* A, double* x, double* y, int n, int it,int n_per_proc, int my_rank);
/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char* name, double*  y, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;

  int mod=0;

  double* local_A=NULL;
  double* local_x=NULL;
  double* local_y=NULL;

  int n, new_n, iters, n_per_proc;
  long seed;
  int my_rank,comm_sz;
  double local_start, local_finish, local_elapsed, elapsed;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  if(my_rank == 0){
    // Obtener las dimensiones
    printf("Ingrese la dimensión n:\n");
    scanf("%d", &n);
    printf("Ingrese el número de iteraciones:\n");
    scanf("%d", &iters);
    printf("Ingrese semilla para el generador de números aleatorios:\n");
    scanf("%ld", &seed);
    
  }
  
  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&seed,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast (&n_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iters,1,MPI_INT,0,MPI_COMM_WORLD);

  n_per_proc = n/comm_sz;
  if(n%comm_sz!=0){
    n_per_proc+=1;
  }
  new_n=n_per_proc*comm_sz;
  A = malloc(sizeof(double) * new_n * new_n);
  srand(seed);
  gen_data(A, new_n*new_n);

  if(my_rank==0){
    x = malloc(sizeof(double) * new_n);
    y = malloc(sizeof(double) * new_n);
    gen_data(x, new_n);
    print_vector("x", x, new_n);
    print_vector("A", A, new_n*new_n);
  }

  local_x=malloc(sizeof(double)*n_per_proc);
  local_y=malloc(sizeof(double)*n_per_proc);
  MPI_Scatter(x, n_per_proc, MPI_DOUBLE, local_x, n_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  char my_rank_str[20];
  sprintf(my_rank_str, "local_x_for:%d", my_rank); 
  print_vector(my_rank_str, local_x, n_per_proc);
  MPI_Barrier(MPI_COMM_WORLD);
  local_start = MPI_Wtime();
  mat_vect_mult(A, local_x, local_y, n, iters,my_rank, n_per_proc, my_rank);
  MPI_Gather(local_y, n_per_proc, MPI_DOUBLE, y, n_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  local_finish = MPI_Wtime();
  local_elapsed = local_finish - local_start;
  
  MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, \
             MPI_MAX, 0, MPI_COMM_WORLD);

  if(my_rank == 0){
    // Solo el proceso 0 imprime el tiempo transcurrido
    printf("Tiempo de ejecución = %5.2f segundos \n", elapsed);
    print_vector("y", y, n);
    free(A);
    free(x);
    free(y);
  }

  MPI_Finalize();
  return 0;
}

void gen_data(double * array, int size){
  int i;
  for (i = 0; i < size; i++)
    array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double* A, double* x, double* y, int n, int it,int n_per_proc, int my_rank){
  int h, i, j;
  for(h = 0; h < it; h++){

    for(i=0; i<n_per_proc; i++){
      y[i] = 0.0;
      for(j = 0; j < n; j++){
         y[i] += A[(my_rank*n_per_proc)+i] * x[i];
      }
    }

    for(i = 0; i < n; i++){
      y[i] = 0.0;
      for(j = 0; j < n; j++){
        y[i] += A[i*n+j] * x[j];
        if(my_rank==0){
          printf("%lf += A[i*n+j]=A[%d]=%lf * x[%d]=%lf=%lf \n",y[i],(i*n+j),A[i*n+j],j,x[j],A[i*n+j] * x[j]);
        }
      }
	      
    }
    // x <= y
    for(i = 0; i < n; i++)
      x[i] = y[i];
  }
}

void print_vector(char* name, double*  y, int m) {
   int i;
   printf("\nVector %s\n", name);
   for (i = 0; i < m; i++){
     printf("%f ", y[i]);
      if(i%16==0){
         printf("\n\n");
      }
   }
      
   printf("\n");
}
