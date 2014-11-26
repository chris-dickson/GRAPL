/*----------------------------------------------------------------------------------*/
/*	globals.h																		*/
/*		Header file for global variables											*/
/*----------------------------------------------------------------------------------*/

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include "binheap.h"
#include "helper.h"


/* Binary heap, and its weight array.   We fill the array position A[i] with the shortest path estimate of vertex i.  
(either 0 for start, or infinity for others),   We then call binheap_build(...)										*/
Heap	*H;
double	*A_weights;

void	init_global_heap(int gw, int gh);	/* allocate space for heap and A_weights (to be 2*gw*gh) */
void	destroy_global_heap();				/* free them */





typedef struct SPElement{
	CoordData 	c;
	double	  	d;
	int		parent;
}SPElement;

int			**xy2term;							/* 2d array that maps terminals (x,y) position to an internal integer representation */
int			*computed_shortest_path;			/* uses internal int representation of a terminal to find if we have the shortest path estimates from that vertex */
SPElement	**shortest_paths;					


/* debug */
int			useless_calls;						/* the number of wasted calls to dijkstra */


void	init_shortest_path_globals(int gw, int gh, int no_unique_terminals);
void	reset_shortest_path(int no_unique_terminals);
int		have_shortest_path(int x, int y);
void	found_shortest_path(int x, int y);		
void	init_xy2term(int gw, int gh, CoordData* unique_terminals, int no_unique_terminals);
void	destroy_shortest_path_globals(int gw, int gh, int no_unique_terminals);

#endif
