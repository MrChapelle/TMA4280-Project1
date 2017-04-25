/* Implementation of the solution to the Poisson problem on [0,1]*[0,1] in an hybrid implementation
 * with a right hand side f . 
 * 
 *  @ authors : You Robin - Houlier Nicolas
 *  @ date    : 23/03/2017
 *  @ TMA 4280 : Introduction to Supercomputing
 *  @ Project 2 Spring 2017
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

#define PI 4.*atan(1.)
#define true 1
#define false 0

typedef double Real;
typedef int Bool;

// Function prototypes

Real *mk_1D_array(size_t n, Bool zero);
Real **mk_2D_array(size_t n1, size_t n2, Bool zero);
Real rhs(Real x, Real y);

// Function related to unit tests 

void show_matrix(double **b, int m, int n);
int is_pow2(int n);

// Functions related to parallelisation
void p_transpose(Real **bt, Real **b, Real *send, Real *recv, int size, int m, int n, int rank);

// Function related to verification of the program
Real sol(Real x, Real y);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.

void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("Not enough arguments:\n");
        printf("./poisson <problem size> <number of threads per process>\n");
    }

    // we initialize the number of processes which will be 
    //the second argument of the function
    
    int t = atoi(argv[2]);
    omp_set_num_threads(t);
    
    // then we initialize the MPI procedure

    int rank;
	int size;
	
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Utilisation de MPI_Wtime pour obtenir le temps courant
    
    Real e_time;
    if (rank == 0) {
      e_time = MPI_Wtime();
    }


    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     *
     * - Number of processors must be a power of two. This allows equal amounts
     *   of data to be handled on each processors making the parallel transpose
     *   simpler. An implementation without this constraint is also be possible.
     */
     
     
    int n = atoi(argv[1]);
    int avg= n/size;
    int nn = 4 * n;
    
    // the step of the grid
    
    Real h = 1.0 / n;
    
    
    // UNIT TESTS 
    
    
    // the main process (0) print the size of the problem , the number of processes , the number of threads and the number of rows per process

    if (rank == 0) printf("Grid size: n = %i, Size = %i, Rows per process: avg= %i, Threads: t = %i\n", n, size, avg, t);

    //Check if size is a power of 2
  	if (!is_pow2(n)) {
      if (rank == 0) {
  			printf("Problem size is not a power of 2. Aborting\n");
  		}
  		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_DIMS);
  		return 0;
    }

	// we define our matrices 
	
    Real **b = mk_2D_array(avg, n, false);
    Real *send = mk_1D_array((size_t)(avg*n), false);
    Real *recv = mk_1D_array((size_t)(avg*n), false);
    Real **bt = mk_2D_array(avg, n, false);
	
	// VERIFICATION TEST :
	// this matrix correspond to the analytic solution
	
    Real **a_sol = mk_2D_array(avg,n, false);

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */

    Real u_max = 0.0, err_max = 0.0;

    Real *x_grid = mk_1D_array(avg, false);
    Real *y_grid = mk_1D_array(n, false);


    // Modification we talked about, z is now 2D not to use
    // previous values for each process
    
    Real **z = mk_2D_array(t,nn,false);


   Real *diag = mk_1D_array(n, false);

    #pragma omp parallel
    {

      #pragma omp for
      for (size_t i = 0; i < avg; i++) {
          x_grid[i] = (i + 1 + avg*rank) * h ;
      }

      #pragma omp for
      for (size_t i = 0; i < n+1; i++) {
          y_grid[i] = (i+1) * h ;
      }

      #pragma omp for
      for (size_t i = 0; i < n; i++) {
          diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));

      }

      #pragma omp for collapse(2)
      for (size_t i = 0; i < avg; i++) {
          for (size_t j = 0; j < n; j++) {
              a_sol[i][j] = sol(x_grid[i], y_grid[j]);
              // h²*f
              b[i][j] = h * h * rhs(x_grid[i], y_grid[j]);
          }
        }

      /*
       * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
       * Instead of using two matrix-matrix products the Discrete Sine Transform
       * (DST) is used.
       * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
       * The array zz is used as storage for DST coefficients and internally for
       * FFT coefficients in fst_ and fstinv_.
       * In functions fst_ and fst_inv_ coefficients are written back to the input
       * array (first argument) so that the initial values are overwritten.
       */

      #pragma omp for
      for (size_t i = 0; i < avg; i++) {
		  
		  // what we talked about : z[]
		  
          fst_(b[i], &n, z[omp_get_thread_num()], &nn);
      }


      #pragma omp master
      p_transpose(bt, b, send, recv, size, avg, n, rank);
      
      #pragma omp barrier
      // identifies a synchronisation point at which point each treads need to wait until
      // the other threads have finished their work

      #pragma omp for
      for (size_t i = 0; i < avg; i++) {
		  // problem we talked about
          fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
      }

      /*
       * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
       */

      #pragma omp for collapse(2)
      for (size_t i = 0; i < avg; i++) {
          for (size_t j = 0; j < n; j++) {
              bt[i][j] = bt[i][j] / (diag[i+rank*avg] + diag[j]);
          }
      }

      /*
       * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
       */

      #pragma omp for
      for (size_t i = 0; i < avg; i++) {
		  // same problem
          fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
      }


      #pragma omp master
      p_transpose(b, bt, send, recv, size, avg, n, rank);
      // identifies a synchronisation point at which point each treads need to wait until
      // the other threads have finished their work
      #pragma omp barrier


      #pragma omp for
      for (size_t i = 0; i < avg; i++) {
          fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
      }

      /*
       * Compute maximal value and maximal absolute error of solution for
       * convergence analysis.
       */

      Real err = 0.0;

      //Last value is not considered since it is on the boundary

      #pragma omp for collapse(2)
        for (size_t i = 0; i < avg; i++) {
          for (size_t j = 0; j < n-1; j++) {
			  // fabs <-> absolute value
              err = fabs(b[i][j] - a_sol[i][j]);
              // the formula you used in your serial program
              err_max = err_max > err ? err_max : err;
              u_max = u_max > b[i][j] ? u_max : b[i][j];
            }
        }
    }

    Real u_max_glob = 0;
    Real err_max_glob = 0;
    
    MPI_Reduce(&u_max, &u_max_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&err_max, &err_max_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
	  // we calculate the final time of the program
      e_time = MPI_Wtime() - e_time ;
      // we print the results
      printf("u_max = %e\n", u_max_glob);
      printf("Time elapsed: %f\n", e_time);
      printf("Max absolute error: %.8f\n\n", err_max_glob);
    }
	
	// we assign the third process to print the final matrix
	// VERIFICATION TEST
    if (rank == 3) show_matrix(b, avg, n-1);

	// We close the MPI routine
    MPI_Finalize();
    
    return 0;
}


/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

Real rhs(Real x, Real y) {
    //return 2 * (y - y*y + x - x*x);
    return 5*PI*PI*sin(PI*x)*sin(2*PI*y); // RHS of known solution
    //return exp(x)*sin(2*PI*x)*sin(2*PI*y);
    //return 1;
}

/*
* This function is used for convergence test of the numerical solution
*/
Real sol(Real x, Real y) {
  return sin(PI*x)*sin(2*PI*y);
}
/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

Real *mk_1D_array(size_t n, Bool zero)
{
    if (zero) {
        return (Real *)calloc(n, sizeof(Real));
    }
    return (Real *)malloc(n * sizeof(Real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

Real **mk_2D_array(size_t n1, size_t n2, Bool zero)
{
    // 1
    Real **ret = (Real **)malloc(n1 * sizeof(Real *));

    // 2
    if (zero) {
        ret[0] = (Real *)calloc(n1 * n2, sizeof(Real));
    }
    else {
        ret[0] = (Real *)malloc(n1 * n2 * sizeof(Real));
    }

    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

void show_matrix(Real **b, int m, int n) {
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
         printf("%f ", b[i][j]);
      }
      printf("\n");
   }
}

/*
* Method that transpose a matrix in parallel. The values in the matrices of
* each process is packed into a 1D array according to the provided figure
* (Figure 0.1) in the problem description. When recived the values from the 1D
* array is then unwrapped and origanized back into matrices.
*/

void p_transpose(Real **bt, Real **b, Real *send, Real *recv, int size, int m, int n, int rank) {
  size_t i, j;
  //#pragma omp parallel for collapse(2)
  for (i=0; i<(size_t)m; i++) {
    for (j=0; j<(size_t)n; j++) {
      send[m*i + (j/m)*(m*m) + j%m] = b[i][j];
    }
  }

  MPI_Alltoall(&send[0], m*m, MPI_DOUBLE, &recv[0], m*m, MPI_DOUBLE, MPI_COMM_WORLD);

  int cnt = 0;
  for (j=0; j<(size_t)n; j++) {
    for (i=0; i<(size_t)m; i++) {
      bt[i][j] = recv[cnt];
      cnt++;
    }
  }
}

int is_pow2(int n) {
  Real l = log(n)/log(2);
  return floor(l) == l;
}
