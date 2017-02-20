/* Authors : You Robin , Houlier Nicolas 
 * Date 08.02.2017
 * This program is the maths test for the first question
 * mach method */

# include <stdio.h>
# include <stdlib.h>

# include <math.h>
# include "question1.h"


double vtest(int k)
{
    FILE *fichier = NULL;
    
    fichier = fopen("vtest.txt","r+");
    
           
    int n = 2;
    double result = 0;
    int i;
    
    for (i = 1;i< (k+1);i++)
    {
		result = fabs(M_PI - question1v2(n));
		printf("La valeur de pi - pin pour n = 2**%d vaut : %.30e \n", i, result);
		n *= 2;
		if (fichier != NULL)
		{
			fprintf(fichier,"pi - pin pour n = 2**%d vaut : \n %.30e \n",i,result);
			
		}
		else
		{
			printf("there is a problem");
		}
	}
	
	fclose(fichier);
	
	return 0;
	
}



