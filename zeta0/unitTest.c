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
	int n = 3;
	double expected_pi = M_PI;
	double piZeta = zeta0sum(n);
	double error = expected_pi - piZeta;
    printf("the error is %f\n", error);
    return 0;
}
