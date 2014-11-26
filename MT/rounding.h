/*--------------------------------------------------------------------------------------*/
/* rounding.h																			*/
/*	This is a header for randomized rounding.  These functions are no longer			*/
/* used in this program.   We call a separate program to round the solutions.			*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/


#ifndef ROUNDING_H_
#define ROUNDING_H_

void display_int(int *x, int K, int l);
void display_double(double *x, int K, int l);
void randomized_round(double *x, int K, int l, int *x_bar);

#endif
