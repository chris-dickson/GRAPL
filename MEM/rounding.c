/*--------------------------------------------------------------------------------------*/
/* rounding.c																			*/
/*	This is an implementation of randomized rounding.									*/
/*--------------------------------------------------------------------------------------*/
#ifndef ROUNDING_C_
#define ROUNDING_C_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#include "rounding.h"

/*--------------------------------------------------------------*/
/*	Purpose: displays x array for all iterations 				*/
/*  Returns: nothing											*/
/*  args   :													*/
/*		(int*) x : our indicator variable						*/
/*		(int) K : number of nets 								*/
/*		(int) l : number of iterations							*/
/*  Other Notes:												*/
/*		We use this to display the integer (rounded) solution	*/
/*--------------------------------------------------------------*/
void display_int(int *x, int K, int l){
	int k;
	
	for (k = 1; k <= K; k++) {
		int i,count;
		
		for (count = 0,i = k; count <= l; i+=K,count++) 
			printf("x[%2d] = %d\t",i,x[i]);
		printf("\n");
	}
	
}


/*--------------------------------------------------------------*/
/*	Purpose: displays x array for a given iteration				*/
/*  Returns: nothing											*/
/*  args   :													*/
/*		(int*) x : our indicator variable						*/
/*		(int) K : number of nets 								*/
/*		(int) l : number of iterations							*/
/*  Other Notes:												*/
/*		We use this to display the fractional solution			*/
/*--------------------------------------------------------------*/
void display_double(double *x, int K, int l){
	int k;
	
	for (k = 1; k <= K; k++) {
		int i,count;
		
		for (count = 0,i = k; count <= l; i+=K,count++) 
			printf("x[%2d] = %f\t",i,x[i]);
		printf("\n");
	}
}


/*--------------------------------------------------------------*/
/*	Purpose: perform a randomized round on x					*/
/*  Returns: (int)* x_bar										*/
/*  args   :													*/
/*		(double*) x : our indicator variable					*/
/*		(int) K : number of nets 								*/
/*		(int) l : number of iterations							*/
/*		(int*) x_bar : our rounded, integer solution			*/
/*--------------------------------------------------------------*/
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
#endif
