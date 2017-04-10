/**
 *  Implementation of the solution to the Poisson problem on [0,1]*[0,1] in a serial way
 *  with a right hand side f . 
 * 
 *  authors : You Robin - Houlier Nicolas
 *  date    : 23/03/2017
 *  TMA 4280 : Introduction to Supercomputing
 *  Project 2 Spring 2017
 *  this code is inspired of the poisson.c program
 * 
 */
 
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>


# define true 1 
# define false 0 

typedef double real;
typedef int bool;

// Functions prototypes
real *make_1D_array (size_t n, bool zero);
real **make_2D_array (int n1, int n2);
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);



int main(int argc, char **argv)
{
	
	if (argc < 2) {
		printf("Problem because n must be a power of 2 \n");
	}
	// Real time= WallTime();
	
	// The number of grid points in each direction is n+1
	// The number of degrees of freedom in each directions is n-1
	
	int n = atoi(argv[1]);
	int m = n - 1;
	int nn = 4*n;
	real h = 1.0/n;
	
	
	// we define our grid of points 
	
	real *grid = make_1D_array(n+1, false);
	
	// we affect the values to our grid of points
	
	for (size_t i = 0; i < n+1; i++) {
		grid[i] = i * h;
	}
	
	// T = Q*LAMBDA*Qt because T is symetrical so can be diagonalised
	// LAMBDA is the matrix of eigenvalue
	// Q is the matrix of eigenvectors [columns]
	
	// we create the diagonal of the eigenvalue of matrix T :
	// T in R(n-1*n-1)
	
	real *diag = make_1D_array(m, false);
	
	for (size_t i = 0; i < m ; i++) {
		diag[i] =  2.0 * (1.0 - cos((i+1) * 3.14159 / n));
	}
	
	// Now we initialise the right hand side , G
	// We have G = h²* Fmatrix 
	// fij = rhs (grid[i] , grid[j]) because rhs define f(x,y)
	
	real **b = make_2D_array(m,m);
	real **bt = make_2D_array(m,m);
	
	// we define the vector z which will contain the Fourier complex 
	// coefficients on -f +f 
	real *z = make_1D_array(nn, false);
	
	// We actualise the coefficient of b
	
	for (size_t i = 0; i<m; i++) {
		for (size_t j=0; j<m;j++) {
			b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
		}
	}
	
	// Now we will do the first step of the resolution
	// We will calculate : GtildeT = QtGQ ie Btilde^T = S^-1 * (S * B)^T
	
	for (size_t i = 0 ; i < m ; i++) {
		fst_(b[i], &n, z, &nn);
	}
	
	// now we transpose SG ie SB :
	
	transpose(bt,b,m);
	
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
	
	for (size_t i=0; i<m;i++) {
		for (size_t j=0; j<m;j++) {
			u_max = u_max > b[i][j] ? u_max : b[i][j];
		}
	}
	
	printf("u_max = %e\n", u_max);
	
	
	return 0;	
	
}


// Other Functions 

real *make_1D_array(size_t n, bool zero)
{
	if (zero) {
		return (real *)calloc(n, sizeof(real));
	}
	return (real *)malloc(n * sizeof(real));
}

real **make_2D_array(int n1, int n2)
{
	real **ret;
	ret = (real **)malloc(n1   *sizeof(real *));
	ret[0] = (real  *)malloc(n1*n2*sizeof(real));
	for (int i=1; i < n1; i++) {
		ret[i] = ret[i-1] + n2;
	}
	int n = n1*n2;
	memset(ret[0],0,n*sizeof(real));
	return ret;
}


void transpose(real **bt, real **b, size_t m)
{
	// the transpose of b is set into bt
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++){
			bt[i][j] = b[j][i];
		}
	}
}

real rhs(real x, real y) {
	return 2 * (y - y*y + x - x*x);
}




