/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * A short test to compute the error between real pi and
 * the pi computed by the main.c programme using zeta method
 * for n=2^k */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "question1.h"


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

