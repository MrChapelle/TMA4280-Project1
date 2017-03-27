/* Implementation of the solution to the Poisson problem on [0,1]*[0,1] in a serial way
 * with a right hand side f . 
 * 
 *  @ authors : You Robin - Houlier Nicolas
 *  @ date    : 23/03/2017
 *  @ TMA 4280 : Introduction to Supercomputing
 *  @ Project 2 Spring 2017
 *  this code is inspired of the poisson.c program
 * 
 */
 
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <memory.h>
# include <sys/time.h>
# include <stddef.h>

# define pi 4.*atan(1.)
# define true 1 
# define false 0 

typedef double Real;
typedef int Bool;

// Functions prototypes

Real *make_1D_array (size_t n, Bool zero);
Real **make_2D_array (size_t n1, size_t n2, Bool zero);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void transpose(Real **bt, Real **b, size_t m);
Real rhs(Real x, Real y);
Real WallTime ();


// Main Function

int main(int argc, char **argv)
{
	if (argc < 2) {
		printf("Problem because n must be a power of 2");
	}
	
	Real time= WallTime();
	
	// The number of grid points in each direction is n+1
	// The number of degrees of freedom in each directions is n-1
	
	int n = atoi(argv[1]);
	int m = n - 1;
	int nn = 4*n;
	Real h = 1.0/n;
	
	// we define our grid of points 
	Real *grid = make_1D_array(n+1, false);
	
	// we affect the values to our grid of points 
	for (size_t i = 0; i < n+1; i++) {
		grid[i] = i * h;
	}
	
	// T = Q*LAMBDA*Qt because T is symetrical so can be diagonalised
	// LAMBDA is the matrix of eigenvalue
	// Q is the matrix of eigenvectors [columns]
	
	// we create the diagonal of the eigenvalue of matrix T :
	// T in R(n-1*n-1)
	
	Real *diag = make_1D_array(m, false);
	
	for (size_t i = 0; i < m ; i++) {
		diag[i] =  2.0 * (1.0 - cos((i+1) * pi / n));
	}
	
	// Now we initialise the right hand side , G
	// We have G = h²* Fmatrix 
	// fij = rhs (grid[i] , grid[j]) because rhs define f(x,y)
	
	Real **b = make_2D_array(m,m,false);
	Real **bt = make_2D_array(m,m,false);
	
	// we define the vector z which will contain the Fourier complex 
	// coefficients on -f +f 
	Real *z = make_1D_array(nn, false);
	
	// We actualise the coefficient of b
	
	for (size_t i = 0; i<m; i++) {
		for (size_t j=0; j<m;j++) {
			b[i][j] = h * h * rhs(grid[i], grid[j]);
		}
	}
	
	// Now we will do the first step of the resolution
	// We will calculate : GtildeT = QtGQ ie Btilde^T = S^-1 * (S * B)^T
	
	for (size_t i = 0 ; i < m ; i++) {
		fst_(b[i], &n, z, &nn);
	}
	
	// now we transpose SG ie SB :
	
	transpose(bt , b , m);
	
	// now we calculate S-1 * bt
	
	for (size_t i = 0; i < m; i++) {
		fstinv_(bt[i], &n , z , &nn);
	}
	
	// the first step is done in O(n²logn)
	
	// now we begin the second step 
	// we want to solve Lambda * Xtilde = Btilde
	// ie : Utilde(ij) = Btilde(ij)/(lj + li)
	
	for (size_t i=0; i<m; i ++) {
		for (size_t j=0; j<m; j++) {
			bt[i][j] = bt[i][j] / ( diag[i] + diag[j]);
		}
	}
	
	// the second step is done in O(n²)
	
	// now we begin the third step
	// calculate U in O(n²logn) ie X= S^-1 * ( S * Xtilde)^T
	
	for (size_t i = 0; i < m ; i++) {
		fst_(bt[i], &n, z, &nn);
	}
	
	transpose(b,bt,m);
	
	for (size_t i = 0; i < m; i++) {
		fstinv_(b[i], &n, z, &nn);
	}
	
	// the third step is done
	
	// we want to return the maximum of our solution U 
	
	double u_max = 0.0;
	
	for (size_t j=0; j<m;j++) {
		for (size_t i=0; i<m;i++) {
			u_max = u_max > b[j][i] ? u_max : b[j][i];
		}
	}
	
	printf("u_max = %e \n", u_max);
	printf("times %f \n ", WallTime()-time);
	
	return 0;	
	
}


// Other Functions 

Real *make_1D_array(size_t n, Bool zero)
{
	if (zero) {
		return (Real *)calloc(n, sizeof(Real));
	}
	return (Real *)malloc(n * sizeof(Real));
}

Real **make_2D_array(size_t n1, size_t n2, Bool zero)
{
	Real **ret = (Real **)malloc(n1 * sizeof(Real *));
	
	if (zero) {
		ret[0] = (Real *)calloc(n1 * n2, sizeof(Real));
	}
	else {
		ret[0] = (Real *)malloc(n1 * n2 * sizeof(Real));
	}
	for (size_t i = 1; 1 < n1 ; i++) {
		ret[i] = ret[i-1] + n2;
	}
	return ret;
}


void transpose(Real **bt, Real **b, size_t m)
{
	// the transpose of b is set into bt
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++){
			bt[i][j] = b[j][i];
		}
	}
}

Real rhs(Real x, Real y)
{
	return 2 * (y - y*y + x - x*x);
}

Real WallTime ()
{
	# ifdef HAVE_MPI
		// we loaded open mpi so we can use the mpi function time
		return MPI_Wtime();
	# endif
	# ifdef HAVE_OPENMP
		// we only loaded openmp so we can use the openmp function time
		return omp_get_Wtime();
	#endif
		// if it is a serial code we take the current time of the day ..
		struct timeval tmpTime;
		gettimeofday(&tmpTime,NULL);
		return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
}



//*****************************END********************************//

	
