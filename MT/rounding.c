/*--------------------------------------------------------------------------------------*/
/* rounding.c																			*/
/*	This is an implementation of randomized rounding.  These functions are no longer    */
/* used in this program.   We call a separate program to round the solutions.			*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#include "rounding.h"

void display_int(int *x, int K, int l){
	int k;
	
	for (k = 1; k <= K; k++) {
		int i,count;
		
		for (count = 0,i = k; count <= l; i+=K,count++) 
			printf("x[%2d] = %d\t",i,x[i]);
		printf("\n");
	}
	
}

void display_double(double *x, int K, int l){
	int k;
	
	for (k = 1; k <= K; k++) {
		int i,count;
		
		for (count = 0,i = k; count <= l; i+=K,count++) 
			printf("x[%2d] = %f\t",i,x[i]);
		printf("\n");
	}
}

void randomized_round(double *x, int K, int l, int *x_bar) {
	int k;
	float rand_num;

	rand_num = (float) rand() / RAND_MAX;
	
	/* for each net*/
	for (k = 1; k <= K; k++) {
		int		i,count;
		double	range;
		
		rand_num = (float) rand() / RAND_MAX;
		range = rand_num;
		
		for (count = 0,i = k; count <= l; i+=K,count++) {
			if ( range - x[i] < 0 ) {
				range = DBL_MAX;
				x_bar[i] = 1;
			}
			else {
				range = range-x[i];
				x_bar[i] = 0;
			}
		}
	}
}		
