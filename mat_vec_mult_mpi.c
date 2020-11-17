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
void mat_vect_mult(double* A, double* x, double* y, int n, int it);
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

  int n, iters, n_per_proc;
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
    n_per_proc = n/comm_sz;
  }
  
  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&seed,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  A = malloc(sizeof(double) * n * n);
  srand(seed);
  gen_data(A, n*n);
  x = malloc(sizeof(double) * (n_per_proc*comm_sz));
  if(my_rank==0){
    y = malloc(sizeof(double) * (n_per_proc*comm_sz));
    gen_data(x, (n_per_proc*comm_sz));
  }
  MPI_Bcast(&iters,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast (&n_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

  local_x=malloc(sizeof(double)*n_per_proc);
  local_y=malloc(sizeof(double)*n_per_proc);

  int AA[16]={1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
  int CC[4];

  MPI_Scatter(AA, 4, MPI_INT, CC, 4, MPI_INT, 0, MPI_COMM_WORLD);

  printf("sale de scater  from process=%d np=%d\n",my_rank,CC[0]);
  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank==0){
    printf("Entra al segundo"); 
  }
  MPI_Scatter(x, n_per_proc, MPI_DOUBLE, local_x, n_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  printf("sale de scater  from process=%d np=%d\n",my_rank,n_per_proc);
  
  
  //Nos aseguramos que todos los procesos inicien al "mismo" tiempo
  MPI_Barrier(MPI_COMM_WORLD);
  printf("continuaaaa p=%d",my_rank);
  local_start = MPI_Wtime();
  printf("from process=%d local_A[0]=%lf\n",my_rank,local_A[0]);
  mat_vect_mult(A, local_x, local_y, n_per_proc, iters);
  MPI_Gather(local_y, n_per_proc, MPI_DOUBLE, y, n_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  local_finish = MPI_Wtime();
  // Cada proceso toma un tiempo local
  local_elapsed = local_finish - local_start;
  // Tomamos el tiempo del proceso más lento
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, \
             MPI_MAX, 0, MPI_COMM_WORLD);

  if(my_rank == 0){
    // Solo el proceso 0 imprime el tiempo transcurrido
    printf("Tiempo de ejecución = %5.2f segundos \n", elapsed);
    //print_vector("y", y, n);
    //free(A);
    //free(x);
    //free(y);
  }

  MPI_Finalize();
  return 0;
}

void gen_data(double * array, int size){
  int i;
  for (i = 0; i < size; i++)
    array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double* A, double* x, double* y, int n, int it){
  int h, i, j;
  for(h = 0; h < it; h++){
    for(i = 0; i < n; i++){
      y[i] = 0.0;
      for(j = 0; j < n; j++)
	      y[i] += A[i*n+j] * x[j];
    }
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
