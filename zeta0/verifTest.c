/*
main.c
@author : Nicolas Houlier 
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "question1.h"

/* A unit test for question1 */

int main(int argc, char *argv[])
{
	double expected_pi = M_PI;
	double error[24] = {0};
	int i = 1;
	for (i = 1 ; i < 25 ; i++)
    {
	double piZeta = zeta0sum(pow(2,i));
	error[i] = expected_pi - piZeta;
    printf("the error is %f\n", error[i]);
    }
    return 0;
}

