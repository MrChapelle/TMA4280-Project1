/* @Author : You Robin & Houlier Nicolas
 * @Date : 19/02/17 
 * A short test to compute the error between real pi and
 * the pi computed by the main.c programme using zeta method*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "question1.h"

int main(int argc, char *argv[])
{
	int n = 3;
	double expected_pi = M_PI;
	double piZeta = zeta0sum(n);
	double error = expected_pi - piZeta;
    printf("the error is %f\n", error);
    return 0;
}
