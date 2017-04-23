/* Implementation of the solution to the Poisson problem on [0,1]*[0,1] in an hybrid implementation
 * with a right hand side f . 
 * 
 *  @ authors : You Robin - Houlier Nicolas
 *  @ date    : 23/03/2017
 *  @ TMA 4280 : Introduction to Supercomputing
 *  @ Project 2 Spring 2017
 * 
 */

# include <mpi.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <memory.h>
# include <stddef.h>
# include <omp.h>


typedef double Real;
typedef int Bool;

# define pi 4.*atan(1.)
# define true 1 
# define false 0 


// functions prototypes
Real *make_1D_array (int n, Bool zero);
Real **make_2D_array (int n1, int n2);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void transpose(Real **Mat, int size, int m, int n, int avg_num , int last_avg_num);
Real rhs(Real x, Real y);

// we are working with mpi so the function walltime is given with mpi_walltime !


// main : 

int main(int argc, char **argv)
{
	
	if (argc < 2) {
		printf("Problem because n must be a power of 2");
	}
	Real **Mat;
	int rank;
	int size;
	int n = atoi(argv[1]);
	
	
	
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	
	/* the total number of grid points in each direction is (n+1)
	 * the total number of degree of freedom in each direction is (n-1)
	*/
	
	int m = n - 1;
	int nn = 4*n;
	Real h = 1.0/n;
	printf("%d \n", size);
	
	int avg_glob = floor(m/size); // nombre moyen de rows par process
	int avg_num = avg_glob * avg_glob; // average number of elements per process
	int last = avg_glob; // in case it's the last process
	int recvcounts = m - (size-1) * avg_glob; 
	int last_avg_num = avg_glob * recvcounts;
	
	if(rank+1 == size) {
    last   = recvcounts;
    avg_num  = last_avg_num;
    last_avg_num = recvcounts * recvcounts;
	}
	
	Real *diag = make_1D_array(m , false);
	
	Real *z = make_1D_array (nn, false);
	
	Mat = make_2D_array (last,m);
	
	Real time = MPI_Wtime();
	
	#pragma omp parallel for schedule(static)
	for (int i=0; i < m; i++) {
		diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
	}
	
	#pragma omp parallel for schedule(static)
	for (int j=0; j < last; j++) {
		for (int i=0; i < m; i++) {
		//        h^2 * f(x,y)
		//Mat[j][i] = h*h*5*pi*pi*sin(pi*i*h)*sin(2*pi*(j + rank*avg_glob)*h);
		Mat[j][i] = h*h*rhs((i+1)*h,(j+1)*h);
		}
	}
	
	#pragma omp parallel for schedule(static)
	for (int j=0; j < last; j++) {
		fst_(Mat[j], &n, z, &nn);
	}

	transpose(Mat, size,last, m, avg_num, last_avg_num);
	
	#pragma omp parallel for schedule(static)
	for (int i=0; i < last; i++) {
		fstinv_(Mat[i], &n, z, &nn);
	}  

	for (int j=0; j < last; j++) {
		for (int i=0; i < m; i++) {
		Mat[j][i] = Mat[j][i]/(diag[i]+diag[j + rank*avg_glob]);
		}
	}
	
	#pragma omp parallel for schedule(static)
	for (int i=0; i < last; i++) {
		fst_(Mat[i], &n, z, &nn);
	}

	transpose(Mat, size, last, m, avg_num, last_avg_num);
	
	#pragma omp parallel for schedule(static)
	for (int j=0; j < last; j++) {
		fstinv_(Mat[j], &n, z, &nn);
	}

	Real umax = 0.0;
	Real emax = 0.0;
	
	for (int j=0; j < last; j++) {
		for (int i=0; i < m; i++) {
		// error =  abs( numerical u(x,y) - exact u(x,y) )
		Real error = fabs(Mat[j][i] - sin(pi*i*h)*sin(2*pi*(j + rank*avg_num)*h));
		if (Mat[j][i] > umax) umax = Mat[j][i];
		if (error > emax) emax = error;
		}
	}
	Real globalumax;
	Real globalemax;

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


// other functions :

void transpose (Real **Mat, int size, int m, int n, int avg_num, int last_avg_num)
{
	// we will follow the hint and use the function MPI_Alltoallv
	// we define our values :
	
	// starting address of send buffer (choice) 
	Real *sendbuf = make_1D_array(n*m, false);
	
	//address of receive buffer (choice) 
	Real *recvbuf = make_1D_array(n*m , false);
	
	// integer array equal to the group size specifying the
	// maximum number of elements that can be received from each processor 
	int recvcounts[size];
	
	// integer array equal to the group size specifying the maximum number of elems 
	// that can be received from each processor 
	int rdispls[size];
	
	// integer array equal to the group size specifying the number of elements to send to each processor
	int sendcounts[size];
	
	// integer array (of length group size). Entry j specifies the displacement (relative to 
	// sendbuf from which to take the outgoing data destined for process j
	int sdispls[size];
	
	// we actualise the values of the parameters :
	
	#pragma omp parallel for schedule(static)
	for (int i=0; i < size; i++) {
		
		sendcounts[i] = avg_num;
		recvcounts[i] = sendcounts[i];
		
		sdispls[i] = avg_num * i ;
		rdispls[i] = sdispls[i];
	}
	
	// we have to be more precise for the last row, the number of rows is not the same
	
	sendcounts[size - 1] = last_avg_num;
	recvcounts[size - 1] = sendcounts[size - 1];
	
	// we spread the information 
	
	#pragma omp parallel for schedule(static)
	for (int i = 0; i<n; i++) {
		for (int j = 0; j < m; j++) {
			sendbuf[i*m + j] = Mat[j][i];
		}
	}
	
	
	// now we apply our MPI_Alltoallv function :
	
	MPI_Alltoallv ( sendbuf , sendcounts , sdispls , MPI_DOUBLE , recvbuf ,
					recvcounts , rdispls , MPI_DOUBLE , MPI_COMM_WORLD);
					
	// now that we have the information 
	// we can actualise the value of Mat
	
	int recvamounts[size];
	
	#pragma omp parallel for schedule(static)
	for (int k = 0; k < size ; k++) {
		recvamounts[k] = recvcounts[k]/m;
	}
	
	// we actualise Mat 
	
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < m ; i++) {
		
		int id = 0;
		int val_int = rdispls[id] + recvamounts[id]*i; // same as line 86
		int val_f = val_int + recvamounts[id] -1;
		
		for (int j = 0; j < n; j++) {
			Mat[i][j] = sendbuf[val_int];
			
			// if it is the last one :
			
			if (val_int == val_f) {
				id ++;
				val_int = rdispls[id] + recvamounts[id]*i; // same as line 110
				val_f = val_int + recvamounts[id] -1; // same as line 111
			}
			else {
				val_int++;
			}
		}
	}
}

Real rhs(Real x, Real y)
{
	return 2 * (y - y*y + x - x*x);
}

Real *make_1D_array(int n, Bool zero)
{
	if (zero) {
		return (Real *)calloc(n, sizeof(Real));
	}
	return (Real *)malloc(n * sizeof(Real));
}

Real **make_2D_array(int n1, int n2)
{
	Real **ret;
	ret = (Real **)malloc(n1   *sizeof(Real *));
	ret[0] = (Real  *)malloc(n1*n2*sizeof(Real));
	for (int i=1; i < n1; i++) {
		ret[i] = ret[i-1] + n2;
	}
	int n = n1*n2;
	memset(ret[0],0,n*sizeof(Real));
	return ret;
}	
	
