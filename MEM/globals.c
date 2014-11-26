/*----------------------------------------------------------------------------------*/
/*	globals.c																		*/
/*		This is the global include file for all other files in the program			*/
/*	We define all global procedures and global variables in the code.				*/
/*----------------------------------------------------------------------------------*/

#ifndef GLOBALS_C_
#define GLOBALS_C_

#include <stdlib.h>
#include "binheap.h"
#include "globals.h"
#include "helper.h"
#include "fatal.h"

/*--------------------------------------------------------------*/
/*	Purpose: initialize a binary heap for comuting shortest		*/
/*			paths   											*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) gw,gh : width and height of MLGraph				*/
/*--------------------------------------------------------------*/
void	init_global_heap(int gw, int gh) {
	H = binheap_init(gw*gh*2);
	A_weights = (double*)malloc(sizeof(double)*((gw*gh*2)+1));
	
	if (A_weights == NULL)
		Error("global A_weights memory allocation error\n");
}

/*--------------------------------------------------------------*/
/*	Purpose: Deallocate the global heap							*/
/*  Returns: nothing										    */
/*  Args   : none												*/
/*		                                 					    */
/*--------------------------------------------------------------*/
void	destroy_global_heap() {
	binheap_destroy(H);
	free(A_weights);
}

/*--------------------------------------------------------------*/
/*	Purpose: Allocates storage for storing shortest paths.		*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) gw, gh : width and height of G                    */
/*		(int) no_unique_terminals : number of vertices in G		*/
/*				that are a terminal in some net of the problem	*/
/*--------------------------------------------------------------*/
void	init_shortest_path_globals(int gw, int gh, int no_unique_terminals) {
	int i;
	
	if ((xy2term = (int**) malloc(sizeof(int*)*gw))==NULL)
		Error("globals xy2term memory allocation error");

	for (i = 0; i < gw; i++)
		if ((xy2term[i] = (int*) malloc(sizeof(int)*gh))==NULL)
			Error("globals xy2term[i] memory allocation error");
	
	if ((computed_shortest_path = (int*) malloc(sizeof(int)*no_unique_terminals))==NULL)
		Error("globals computed_shortest_path memory allocation error");
	
	/* one list for each unique terminals */
	if ((shortest_paths = (SPElement**) malloc(sizeof(SPElement*)*no_unique_terminals))==NULL)
		Error("globals shortest_paths memory allocation error");
	
	for (i = 0; i < no_unique_terminals; i++)
		/* each array has n elements */
		if ((shortest_paths[i] = (SPElement*) malloc(sizeof(SPElement)*(2*gw*gh)))==NULL)
			Error("shortest_paths[i] memory allocation error");
	
	reset_shortest_path(no_unique_terminals);
}

/*--------------------------------------------------------------*/
/*	Purpose: Call when we compute shortest path for (x,y)		*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) x,y : coordinate of terminal						*/
/*--------------------------------------------------------------*/
void	found_shortest_path(int x, int y) {
	computed_shortest_path[ xy2term[x][y] ] = 1;
}

/*--------------------------------------------------------------*/
/*	Purpose: Clear the shortest path data						*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) no_unique_terminals :  same as above              */
/*	Other )Notes:												*/
/*		This procedure marks the paths as uncomputed.  It does	*/
/*	not destroy the old paths.  This needs to be called each	*/
/*	time the weights in the graph change						*/
/*--------------------------------------------------------------*/
void	reset_shortest_path(int no_unique_terminals) {
	int i;
	
	for (i = 0; i < no_unique_terminals; i++)
		computed_shortest_path[i] = 0;
}


/*--------------------------------------------------------------*/
/*	Purpose: Checks if we have shortest path for (x,y)			*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) x,y :  coordinates of terminal                    */
/*--------------------------------------------------------------*/
int		have_shortest_path(int x, int y) {
	return( computed_shortest_path[ xy2term[x][y] ] );
}


/*--------------------------------------------------------------*/
/*	Purpose: Initializes lookup array.  						*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) gw,gh : width and height of G                     */
/*		(CoordData*) unique_terminals : array of unique			*/
/*						terminals								*/
/*		(int) no_unique_terminals : size of above argument		*/
/*--------------------------------------------------------------*/
void	init_xy2term(int gw, int gh, CoordData *unique_terminals, int no_unique_terminals) {
	int	i,j;
	
	/* initialize so that (x,y)'s that aren't terminals return -1 */
	for (i = 0; i < gw; i++) 
		for (j = 0; j < gh; j++) 
			xy2term[i][j] = -1;
	
	/* set the index */
	for (i = 0; i < no_unique_terminals; i++) 
		xy2term[ unique_terminals[i].x ][ unique_terminals[i].y ] = i;
}

/*--------------------------------------------------------------*/
/*	Purpose: Frees all global variables							*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) gw,gh : width and height of G                     */
/*		(int) no_unique_terminals : same as above				*/
/*--------------------------------------------------------------*/
void	destroy_shortest_path_globals(int gw, int gh, int no_unique_terminals) {
	int i;
	
	for (i = 0; i < no_unique_terminals; i++) 
		free(shortest_paths[i]);
	free(shortest_paths);
	
	for (i = 0; i < gw; i++)
		free(xy2term[i]);
	free(xy2term);
	
	free(computed_shortest_path);
}

#endif
