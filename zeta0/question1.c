/*
main.c
@author : Nicolas Houlier 
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "question1.h"

/* A function which compute pi */

double zeta0sum(int n)
{
	// first lets initialize our variables, sum represents the serie Sn
	int count;
	double Sn = 0;
	
	// Lets compute the sum of 1/i^2 from 0 to n
	for (count = 1 ; count < n ; count ++)
	{
		Sn = Sn + 1./(double)(count*count);
	}
	// Since the limit for infite n of Sn = pi^2/6, we compute pi from Sn
    return sqrt(6*Sn);
}    
