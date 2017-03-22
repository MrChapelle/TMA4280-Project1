/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * Programm which calculate the sum of an array using open MPI 
 * The root process acts as a master
 * He sends a portion of the array to each process
 * Master and child calculate their partial sum
 * Children send it to the root process
 * Root process calculate the total sum 
 * Root process print it */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include <mpi.h>
# include <omp.h>


# define max_rows 600000
# define send_data_tag 2001
# define return_data_tag 2002


double array[max_rows];
double array2[max_rows];

double main (int argc , char **argv)
{
	printf("Enter the number of numbers to sum \n");
	    
	
	double sum , partial_sum ;
	
	double result = 0;
	double error = 0;
	
	double x1 = (double)1/5;
	double x2 = (double)1/239;
	
	int num_rows;
	
	/* I query how many number I will sum in order to distribute the work */
		
	
	scanf("%i", &num_rows);
	
			
	MPI_Status status;
	
	int my_id , root_process , ierr , i  , num_procs ,
	    an_id , num_rows_to_receive , avg_rows_per_process , sender , 
	    num_rows_received , start_row , end_row , num_rows_to_send ;
	    
	int num_procs2;
	
	
	
	/* Now each process will execute a separate copy of this program */
	
	ierr = MPI_Init(&argc,&argv);
	root_process = 0 ;
	
	/* Want my current process id and the number of processes started */
	
	ierr = MPI_Comm_rank(MPI_COMM_WORLD , &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD , &num_procs);
	
	
	if (!(num_procs == 1) && !(num_procs == 2) && !(num_procs == 4) &&
	    !(num_procs == 8) && !(num_procs == 16) && !(num_procs == 32) &&
	    !(num_procs == 64) && !(num_procs == 128) && !(num_procs == 256))
	{
		printf("invalid number of processes, not a power of 2 ... \n");
		printf("%i%", num_procs);
		exit(1);
	}
	
	
	
	/* first case : I am the root process : */
	
	if (my_id == root_process)
	{
				
		/* I verify that this number is not too big */
		
		if (num_rows > max_rows)
		{
			printf("Too many numbers \n");
			exit(2);
		}
		
		
		avg_rows_per_process = num_rows / num_procs ;
		
		/* creation of the vector to sum */
		
		for (i = 0 ; i < num_rows ; i++)
		{
			array[i] = (4*(pow(-1,i)*pow(x1,2*i+1))-pow(-1,i)*pow(x2,2*i+1))/(2*i+1);
			
		}			
	
	    /* Now we distribute a portion of the vector to each child process */
	    
	    for(an_id = 1 ; an_id < num_procs ; an_id ++)
	    {
			start_row = an_id * avg_rows_per_process;
			end_row = (an_id + 1) * avg_rows_per_process - 1;
			
			if(an_id == num_procs - 1)
			{
				end_row = num_rows +1 ;
			}
			
			num_rows_to_send = end_row - start_row + 1;
			
			ierr = MPI_Send(&num_rows_to_send, 1 , MPI_INT , an_id , send_data_tag , MPI_COMM_WORLD);
			ierr = MPI_Send(&array[start_row], num_rows_to_send , MPI_DOUBLE , an_id , send_data_tag , MPI_COMM_WORLD);
			
		}
		
		/* then we calculate the sum of the values affected to the root process */
		
		/********* we change this part using open MP **************/
		
		# pragma omp parallel shared (array, num_procs2 , avg_rows_per_process , sum) private(i)
		{
			
			sum = 0;
			num_procs2 = omp_get_num_threads();
		
			# pragma omp for schedule(static)
		
			for (i=0;i<avg_rows_per_process ; i++)
			{
				sum += array[i];
			}
		
			printf("sum %e calculated by the root process \n",sum);
		}
		
		/***********************************************************/
		
		/* then we collect the partial sums */
		
		for (an_id = 1 ; an_id < num_procs ; an_id ++)
		{
			ierr = MPI_Recv(&partial_sum, 1 , MPI_DOUBLE , MPI_ANY_SOURCE , return_data_tag , MPI_COMM_WORLD , &status);
			
			sender = status.MPI_SOURCE;
			printf("Partial sum %e returned from process %i \n", partial_sum , sender );
			
			sum += partial_sum ;
		}
		
				
		result = 4*sum;
		
		printf("Pi is approached by %.30e \n" , result);
		error = fabs(M_PI - result);
		printf("error : %.30e \n", error);	
		
		return 0;
		
		
			
	} 
	
	else
	{
		ierr = MPI_Recv(&num_rows_to_receive , 1 , MPI_INT , root_process , send_data_tag , MPI_COMM_WORLD, &status);
		ierr = MPI_Recv(&array2 , num_rows_to_receive , MPI_DOUBLE , root_process , send_data_tag , MPI_COMM_WORLD , &status);
		
		num_rows_received = num_rows_to_receive ;
		
		partial_sum = 0;
		
		/************ change this part using open MP *****************/
		
		# pragma omp parallel shared (array2, num_procs2 , num_rows_received, sum) private(i, partial_sum)
		{
			num_procs2= omp_get_num_threads();
			
			# pragma omp for schedule(static)	
		
			for(i=0; i<num_rows_received; i++)
			{
				partial_sum += array2[i];
			}
		
		}
		
		/*************************************************************/
		
		/* we send it to the root process */
		
		ierr = MPI_Send(&partial_sum , 1 , MPI_DOUBLE , root_process , return_data_tag , MPI_COMM_WORLD);
	}
	
	ierr = MPI_Finalize();
	
		
		
}
	
	

