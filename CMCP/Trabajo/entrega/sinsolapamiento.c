#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NREPS 10000

/* 
 * Multiplicación de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 */
void matvec(int N,int b,double *A, double *v, double *w)
{
  int i, j, li, ls;

  for (i=0; i<N; i++) {
    w[i] = 0.0;
    li = i-b<0? 0: i-b;  /* limite inferior */
    ls = i+b>N-1? N-1: i+b;  /* limite superior */
    for (j=li; j<=ls; j++) {
      w[i] += A[i*N+j]*v[j];
    }
  }
}

void matvec_parallel(int N,int nLocal, int b,double *A, double *v, double *w)
{
  int i, j, li, ls,size,rank,next,prev,jglobal;
  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;
  if (rank == 0) prev = MPI_PROC_NULL;
  else prev = rank-1;
  if (rank == size-1) next = MPI_PROC_NULL;
  else next = rank+1;
  MPI_Sendrecv(&v[b], b, MPI_DOUBLE, prev, 0, &v[nLocal+b], b, MPI_DOUBLE, next, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  MPI_Sendrecv(&v[nLocal], b, MPI_DOUBLE, next, 0, &v[0], b, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  for (i=0; i<=nLocal-1; i++) {
    w[i] = 0.0;
    jglobal = (rank*nLocal)+i;
    li = jglobal-b<0? 0: jglobal-b;  /* limite inferior */
    ls = jglobal+b>N-1? N-1: jglobal+b;  /* limite superior */
    for (j=li; j<=ls; j++) {
      w[i] += A[i*N+j]*v[j-(rank*nLocal)+b];
    }
  }
  
}

int main(int argc, char **argv) 
{
  int i, j, k, N=20, b=4, size, rank, nLocal,iglobal,jglobal;
  double *A, *v, *w, *aux,t1,t2;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
     if ((N = atoi(argv[1])) < 0) N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de b */
     if ((b = atoi(argv[2])) < 0) b = 1;
  }
  if (b>=N) { /* Valor de b incorrecto */
    printf("Error: ancho de banda excesivo, N=%d, b=%d\n", N, b);
    exit(1);
  }
  nLocal= N/size;
  t1 = MPI_Wtime();
  
  /* Reserva de memoria */
  A = (double*)calloc(nLocal*N,sizeof(double));
  v = (double*)calloc(nLocal+(2*b),sizeof(double));
  w = (double*)calloc( nLocal,sizeof(double));
  aux = (double*)calloc( N,sizeof(double));

  /* Inicializar datos */
  for (i=0; i<nLocal; i++){
    iglobal = (rank*nLocal)+i;

    A[i*N+iglobal] = 2*b;
  } 
  for (i=0; i<nLocal; i++) {
    iglobal = (rank*nLocal)+i;
    for (j=0; j<N; j++) {
      if (iglobal!=j && abs(iglobal-j)<=b) A[i*N+j] = -1.0;
    }
  }

  for (i=b; i<nLocal+b; i++) v[i] = 1.0;

  
  /* Multiplicación de matrices */
  for (k=0; k<NREPS; k++) matvec_parallel(N,nLocal,b,A,v,w);
 
  /* Imprimir solución */
  MPI_Gather(&w[0], nLocal, MPI_DOUBLE, aux, nLocal,MPI_DOUBLE,0,MPI_COMM_WORLD);
  t2 = MPI_Wtime();

   if(rank==0){
    if (N<100) for (i=0; i<N; i++) printf("w[%d] = %g\n", i, aux[i]);
    printf("El tiempo es %g\n", t2-t1);
   }
  

  free(A);
  free(v);
  free(w);
  free(aux);
  MPI_Finalize();

  return 0;
}