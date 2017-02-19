/* Authors : You Robin , Houlier Nicolas 
 * Date 08.02.2017
 * This program is unit test for the first question
 * mach method */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "question1.h"

double utest(int n)
{
	double expected_result  = 3.141621;
	
	double calculated_result = question1v2(n);
		
	if (fabs(expected_result - calculated_result)>0.000001)
	{
			printf("There is a problem somewhere \n");
	}
	else 
	{
		printf("The program seems to be ok for n = %d \n",n);
	}
	
	return 0;
	
}


