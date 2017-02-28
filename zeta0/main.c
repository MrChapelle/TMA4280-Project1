/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * Programm which compute pi with the zeta method
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "question1.h"

int main(int argc, char *argv[])
{
    int inputNumber = 0;
    double pival = 0;
    
    printf("Enter a number... ");
    scanf("%d", &inputNumber);
    
	// Computation of pi from the number input by the user
    pival = zeta0sum(inputNumber);
    printf("\nthe approached value of pi is %f\n", pival);
    
    return 0;
}



