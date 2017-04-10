/* Implementation of the solution to the Poisson problem on [0,1]*[0,1] in a hybrid implementation
 * OPEN MP & MPI
 * with a right hand side f . 
 * 
 *  @ authors : You Robin - Houlier Nicolas
 *  @ date    : 23/03/2017
 *  @ TMA 4280 : Introduction to Supercomputing
 *  @ Project 2 Spring 2017
 * 
 */


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void transpose (Real **A, int m, int n, int size, int bb, int bre);
void fillA (Real **A, Real *V, int *re, int *rd, int m, int n, int size);

int main(int argc, char **argv )
{
  Real *diag, **A;
  Real pi, h, umax, globalumax, emax, globalemax, error, time;
  int n, m, nn, b, re, l, bb, bre, rank, size;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 if( argc < 2 ) {
    printf("need a problem size\n");
	return 0;
  }

  n  = atoi(argv[1]);
  m  = n-1;
  nn = 4*n;

  h    = 1./(Real)n;
  pi   = 4.*atan(1.);

  b   = floor(m/size);
  re  = m - (size-1)*b;
  l   = b;
  bb  = b*b;
  bre = b*re;

  if(rank+1 == size) {
    l   = re;
    bb  = bre;
    bre = re*re;
  }

  diag = createRealArray (m);
  A    = createReal2DArray (l,m);

  time = MPI_Wtime();

  #pragma omp parallel for schedule(static)
  for (int i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }

  #pragma omp parallel for schedule(static)
  for (int j=0; j < l; j++) {
    for (int i=0; i < m; i++) {
      //        h^2 * f(x,y)
      A[j][i] = h*h*5*pi*pi*sin(pi*i*h)*sin(2*pi*(j + rank*b)*h);
    }
  }
  
  #pragma omp parallel for schedule(static)
  for (int j=0; j < l; j++) {
    Real *z = createRealArray (nn);
    fst_(A[j], &n, z, &nn);
  }

  transpose(A, l, m, size, bb, bre);

  #pragma omp parallel for schedule(static)
  for (int i=0; i < l; i++) {
    Real *z = createRealArray (nn);
    fstinv_(A[i], &n, z, &nn);
  }  

  for (int j=0; j < l; j++) {
    for (int i=0; i < m; i++) {
      A[j][i] = A[j][i]/(diag[i]+diag[j + rank*b]);
    }
  }
  
  #pragma omp parallel for schedule(static)
  for (int i=0; i < l; i++) {
    Real *z = createRealArray (nn);
    fst_(A[i], &n, z, &nn);
  }

  transpose(A, l, m, size, bb, bre);

  #pragma omp parallel for schedule(static)
  for (int j=0; j < l; j++) {
    Real *z = createRealArray (nn);
    fstinv_(A[j], &n, z, &nn);
  }

  umax = 0.0;
  emax = 0.0;
  for (int j=0; j < l; j++) {
    for (int i=0; i < m; i++) {
      // error =  abs( numerical u(x,y) - exact u(x,y) )
      error = fabs(A[j][i] - sin(pi*i*h)*sin(2*pi*(j + rank*b)*h));
      if (A[j][i] > umax) umax = A[j][i];
      if (error > emax) emax = error;
    }
  }

  MPI_Reduce (&umax, &globalumax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (&emax, &globalemax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (rank == 0)
  {
    printf("elapsed: %f\n", MPI_Wtime()-time);
    printf ("umax = %e \n",globalumax);
    printf ("emax = %e \n",globalemax);
  }

  MPI_Finalize();
  return 0;
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}

void transpose (Real **A, int m, int n, int size, int bb, int bre)
{
  int se[size], sd[size], re[size], rd[size];
  Real *V = createRealArray (n*m);
  Real *Vt = createRealArray (n*m);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < size; ++i) {
    se[i] = bb;
    sd[i] = bb*i;
    re[i] = bb;
    rd[i] = bb*i;
  }
  se[size-1] = bre;
  re[size-1] = bre;

  #pragma omp parallel for schedule(static)
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      V[j + i*m] = A[j][i];
    }
  }

  MPI_Alltoallv(V, se, sd, MPI_DOUBLE, Vt, re, rd, MPI_DOUBLE, MPI_COMM_WORLD);
  fillA(A, Vt, re, rd, m, n, size);
}

void fillA (Real **A, Real *V, int *re, int *rd, int m, int n, int size)
{
  int rem[size];

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < size; ++i) {
    rem[i] = re[i]/m;
  }

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < m; i++) {
    int r, k, k1;

    r  = 0;
    k  = rd[r] + rem[r]*i;
    k1 = k + rem[r] - 1;

    for (int j = 0; j < n; j++) {
      A[i][j] = V[k];

      if(k == k1) {
        r++;
        k  = rd[r] + rem[r]*i;
        k1 = k + rem[r] - 1;
      } else {
        k++;
      }
    }
  }
}
