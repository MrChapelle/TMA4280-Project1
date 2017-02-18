/* Authors : You Robin , Houlier Nicolas 
 * Date 08.02.2017
 * 
 * This program solves the question 1 of the Project 1 
 * using Method 2 [Machin Formula] */
 
/* import of libraries
 * n will not exceed a maximum max_rows */
 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "question1.h"
# define max_rows 500000
 
double question1()
{
	int n = 0;
	int i = 0;
	
	double sum1 = 0;
	double sum2 = 0;
	
	double array[max_rows];
	double array2[max_rows];
	
	double x1 = (double)1/5;
	double x2 = (double)1/239;
	
	double result;
	
	
	printf("enter the number n in order to calculate the sum : \n");
	scanf("%d",& n);
	
	/* we check that the number is not too big */
	
	if(n>max_rows)
	{
		printf("This number is too big \n");
		
	}
	
	/* we create our arrays */
	 
	else
	{
		for ( i = 0 ; i < n ; i++)
		{
			array[i] = pow(-1,i)*pow(x1,2*i+1)/(2*i+1);
			array2[i] = pow(-1,i)*pow(x2,2*i+1)/(2*i+1);
			
		}
	}
	
	/* we calculate the sum of our arrays */
	
	for (i=0 ; i < n ; i++)
	{
		sum1 += array[i];
		sum2 += array2[i];
		
	}
	
	result = 4*(4*sum1 -sum2);
	
	printf("the approximation of pi is %.15e : \n", result); 
	
	return 0;
	
		
}

double question1v2(int n)
{
	int i = 0;
	
	double sum1 = 0;
	double sum2 = 0;
	
	double array[max_rows];
	double array2[max_rows];
	
	double x1 = (double)1/5;
	double x2 = (double)1/239;
	
	double result;
	
		
	/* we check that the number is not too big */
	
	if(n>max_rows)
	{
		printf("This number is too big \n");
		
	}
	
	/* we create our arrays */
	 
	else
	{
		for ( i = 0 ; i < n ; i++)
		{
			array[i] = pow(-1,i)*pow(x1,2*i+1)/(2*i+1);
			array2[i] = pow(-1,i)*pow(x2,2*i+1)/(2*i+1);
			
		}
	}
	
	/* we calculate the sum of our arrays */
	
	for (i=0 ; i < n ; i++)
	{
		sum1 += array[i];
		sum2 += array2[i];
		
	}
	
	result = 4*(4*sum1 -sum2);
	
	printf("the approximation of pi is %.15e : \n", result); 
	
	return result;
	
		
}




