/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * Programm which calculate the sum of an array using open MP
 * The root process acts as a master
 * He sends a portion of the array to each process
 * Master and child calculate their partial sum
 * Root process calculate the total sum 
 * Root process print it */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include <omp.h>

# define max_rows 600000


double array[max_rows];
double error;


double main (int argc , char **argv)
{
	printf("Enter the number of numbers to sum \n");
	
	double x1 = (double)1/5;
	double x2 = (double)1/239;
	double sum;
	double result;
	int num_rows;
	
	
	scanf("%i", &num_rows);
	
		
	int i , my_id , num_procs , start_row , end_row , avg_rows_per_process ;
	
	
	sum = 0;
	
	if (num_rows > max_rows)
		{
			printf("Too many numbers \n");
			exit(2);
		}
	
	
	
	for (i = 0 ; i < num_rows ; i++)
		{
			array[i] = 1./(double)((i+1)*(i+1));
			
		}	
		
	# pragma omp parallel shared(array, num_rows, num_procs,sum) private( my_id, i)
	{
		num_procs = omp_get_num_threads();
				
						
		if (!(num_procs == 1) && !(num_procs == 2) && !(num_procs == 4) &&
			!(num_procs == 8) && !(num_procs == 16) && !(num_procs == 32) &&
			!(num_procs == 64) && !(num_procs == 128) && !(num_procs == 256))
		{
			printf("invalid number of processes, not a power of 2 ... \n");
			printf("%i%", num_procs);
			exit(1);
		}
		
		my_id = omp_get_thread_num();
		
		# pragma omp for schedule(static)
		
		for(i=0;i<num_rows;i++)
		{
			sum += array[i];
		}
		
	}
	
	result = sqrt(6*sum);
		
	printf("Pi is approached by %.30e \n" , result);
	error = fabs(M_PI - result);
	printf("error : %.30e \n", error);		
	
}


