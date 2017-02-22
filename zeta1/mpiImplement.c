/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * Programm which calculate the sum of an array using open MPI 
 * The root process acts as a master
 * He sends a portion of the array to each process
 * Master and child calculate their partial sum
 * Children send it to the root process
 * Root process calculate the total sum 
 * Root process print it */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define max_rows 600000
#define send_data_tag 2001
#define return_data_tag 2002


double array[max_rows];
double array2[max_rows];

double main(int argc, char **argv) 
{	
	// We ask the user for a number of terms he wants to sum
	printf("\nEnter a number ... \n");
	
	// And we initialize our variables
	double approachedPi = 0;
	double error = 0;
	double sum, partial_sum;
	int my_id, root_process, ierr, i, num_rows, num_procs, 
		an_id, num_rows_to_receive, avg_rows_per_process, sender, 
		num_rows_received, start_row, end_row, num_rows_to_send;
	
	scanf("%i", &num_rows);

	
	MPI_Status status;
	ierr = MPI_Init(&argc, &argv);
	root_process = 0;

	// find out MY process ID, and how many processes were started

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if (!(num_procs == 1) && !(num_procs == 2) && !(num_procs == 4) &&
		!(num_procs == 8) && !(num_procs == 16) && !(num_procs == 32) &&
		!(num_procs == 64) && !(num_procs == 128) && !(num_procs == 256))
	{
		printf("invalid number of processes, not a power of 2 ... \n");
		printf("%i%", num_procs);
		exit(1);
	}

	if(my_id == root_process) {

		if(num_rows > max_rows) {
			printf("Too many numbers.\n");
			exit(2);
		}

		avg_rows_per_process = num_rows / num_procs;

		// We put each terms in the array

		for(i = 0; i < num_rows; i++) {
			array[i] = 1./(double)((i+1)*(i+1));
		}

		/* We distribute a portion of the bector to each child process */

		for(an_id = 1; an_id < num_procs; an_id++) {
			
			start_row = an_id*avg_rows_per_process + 1;
			end_row   = (an_id + 1)*avg_rows_per_process;

			if((num_rows - end_row) < avg_rows_per_process)
			{
				end_row = num_rows - 1;
			}
			
			num_rows_to_send = end_row - start_row + 1;

			ierr = MPI_Send( &num_rows_to_send, 1 , MPI_INT,
				  an_id, send_data_tag, MPI_COMM_WORLD);

			ierr = MPI_Send( &array[start_row], num_rows_to_send, MPI_INT,
				  an_id, send_data_tag, MPI_COMM_WORLD);
		}

		// We and calculate the sum of the values in the segment assigned
		 *   * to the root process */

		sum = 0;
		
		for(i = 0; i < avg_rows_per_process + 1; i++) 
		{
			sum += array[i];   
		} 

		//printf("sum %i calculated by root process\n", sum);

		/* and, finally, I collet the partial sums from the slave processes, 
		* print them, and add them to the grand sum, and print it */

		for(an_id = 1; an_id < num_procs; an_id++) {
			
			ierr = MPI_Recv( &partial_sum, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
				  return_data_tag, MPI_COMM_WORLD, &status);

			sender = status.MPI_SOURCE;

			// printf("Partial sum %i returned from process %i\n", partial_sum, sender);
	 
			sum += partial_sum;
		}
		
		approachedPi = sqrt(6*sum);
		printf("The approached value of pi is: %.30e\n", approachedPi);
		error = fabs(M_PI - approachedPi);
		printf("error : %.30e \n", error);
	}

	else {

	/* I must be a slave process, so I must receive my array segment,
		* storing it in a "local" array, array1. */

		ierr = MPI_Recv( &num_rows_to_receive, 1, MPI_INT, 
		   root_process, send_data_tag, MPI_COMM_WORLD, &status);

		ierr = MPI_Recv( &array2, num_rows_to_receive, MPI_INT, 
		   root_process, send_data_tag, MPI_COMM_WORLD, &status);

		num_rows_received = num_rows_to_receive;

		/* Calculate the sum of my portion of the array */

		partial_sum = 0;
		for(i = 0; i < num_rows_received; i++) 
		{
			partial_sum += array2[i];
		}

		/* and finally, send my partial sum to hte root process */

		ierr = MPI_Send( &partial_sum, 1, MPI_DOUBLE, root_process, 
		return_data_tag, MPI_COMM_WORLD);
	}

ierr = MPI_Finalize();
}
