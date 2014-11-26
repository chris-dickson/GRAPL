/*----------------------------------------------------------------------------------*/
/*	globals.MT.c																	*/
/*		This is the global include file for all other files in the program			*/
/*	We define all global procedures and global variables in the code.				*/
/*																					*/
/*		THIS FILE IS FOR MULTI-THREADED CODE										*/
/*----------------------------------------------------------------------------------*/


#ifndef GLOBALS_MT_C_
#define GLOBALS_MT_C_

#include <pthread.h>
#include <stdlib.h>
#include "binheap.h"
#include "globals.MT.h"
#include "helper.h"
#include "fatal.h"
#include "steiner.MT.h"
#include "graphtype.h"

/*--------------------------------------------------------------*/
/*	Purpose: Initializes N heaps where N is the number of		*/
/*		threads to be run										*/
/*  Returns: Heap array H, array of double arrays A				*/
/*  Args   :												    */
/*		(Heap ***) H  : pointer to array of Heap pointers		*/
/*		(double***) A : pointer to array of double arrays		*/
/*		(int) gw,gh	  : width and height of G					*/
/*		(int) no_threads : number of threads					*/
/*--------------------------------------------------------------*/
void	init_threaded_heaps(Heap ***H, double ***A, int gw, int gh, int no_threads) {
	int i;
	
	*H = (Heap**) malloc(sizeof(Heap*)*no_threads);
	*A = (double**) malloc(sizeof(double*)*no_threads);
	
	for (i = 0; i < no_threads; i++) {
		(*H)[i] = binheap_init(gw*gh*2);
		(*A)[i] = (double*)malloc(sizeof(double)*((gw*gh*2)+1));
		if ( ((*A)[i]) == NULL ) 
			Error("global A[i] memory allocation error\n");
	}	
}

/*--------------------------------------------------------------*/
/*	Purpose: Destroys globals heap variables					*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) no_threads : number of threads					*/
/*--------------------------------------------------------------*/
void	destroy_threaded_heaps(int no_threads) {
	int i;
	
	for (i = 0; i < no_threads; i++) {
		binheap_destroy(H[i]);
		free(A_weights[i]);
	}
	
	free(H);
	free(A_weights);
}

/*--------------------------------------------------------------*/
/*	Purpose: Initizlies the variables for storing shortest		*/
/*		paths.  We also initialize the synchronization vars		*/
/*		for these variables.									*/
/*  Returns: nothing										    */
/*  Args   :												    */
/*		(int) gw,gh : width and height of G						*/
/*		(int) no_unique_terminals : number of terminals that	*/
/*			exists in some net of the problem					*/
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
	
	if ((computed_shortest_path_lock = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t)*no_unique_terminals))==NULL)
		Error("globals computed_shortest_path_lock memory allocation error");
		
	if ((shortest_paths_lock = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t)*no_unique_terminals))==NULL)
		Error("globals computed_shortest_path_lock memory allocation error");
	
	for (i = 0; i < no_unique_terminals; i++)  {
		pthread_mutex_init(&(computed_shortest_path_lock[i]),NULL);
		pthread_mutex_init(&(shortest_paths_lock[i]),NULL);
	}
		
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
	
	/* must free the locks here too!!*/
}

/*--------------------------------------------------------------*/
/*	Purpose: Initialize the structure that contains all data to	*/
/*		be passed to the worker thread.   Since worker threads	*/
/*		only take a single argument, we have to package all data*/
/*		neccessary to generate a set of trees					*/
/*  Returns: (thread_args_t**) args								*/
/*  Args   :												    */
/*		(thread_args_t**) args : array of arguments for workers	*/
/*		(int) no_threads : number of threads					*/
/*		(Heap**) H : array of Heaps								*/
/*		(double**) A : array of weights							*/
/*		(Graph*) g : ML Graph									*/
/*		(AdjList****) : Edge2Adj : (x,y,z) -> adjacency list	*/
/*		(int) width,height : width and height of ML Graph		*/
/*		(double*) L : vector of edge lengths (index by edges)	*/
/*		(int) K : number of nets								*/
/*		(int) l : iteration number								*/
/*--------------------------------------------------------------*/
void init_thread_args(thread_args_t** args, int no_threads, Heap **H, double **A, Graph **g, AdjList ****Edge2Adj, int width, int height, double *L,  int K,int l) {
	int i;
	
	(*args) = (thread_args_t*) malloc(sizeof(thread_args_t)*no_threads);
	for (i = 0; i < no_threads; i++) {
		thread_args_t *a;
		
		a = &((*args)[i]);
		
		/* set arguments */
		a->H = H[i];
		a->A = A[i];
		a->g = g[i];
		a->Edge2Adj = Edge2Adj;
		a->width = width;
		a->height = height;
		a->L = L;
		a->l = l;
		a->K = K;
		
		/* we now have to split T into equal piles.*/
		a->lb = i*(K/no_threads)+1;	
		a->ub = (i+1)*(K/no_threads);
		
		/* the last thread will pick up the remainder */
		if (i == no_threads-1)
			a->ub += K%no_threads;
		
		a->no_jobs = a->ub - a->lb + 1;
	}

}

/*--------------------------------------------------------------*/
/*	Purpose: worker thread.   This function computes a range of	*/
/*		trees specified by its input argument					*/
/*  Returns: nothing  (modifies T, and shortest paths)			*/
/*  Args   :												    */
/*		(thread_args_t*) args : packaged arguments				*/
/*--------------------------------------------------------------*/
void *runner(void *args) {
	Heap	*H;
	double	*A;
	int		lb,ub,no_jobs;
	Graph	*g;
	AdjList	****Edge2Adj;
	int		width,height;
	double	*L;
	int		l;
	int		K;
	
	int		k;
	
	H = ((thread_args_t*) args)->H;
	A = ((thread_args_t*) args)->A;
	lb = ((thread_args_t*) args)->lb;
	ub = ((thread_args_t*) args)->ub;
	no_jobs = ((thread_args_t*) args)->no_jobs;
	g = ((thread_args_t*) args)->g;
	Edge2Adj = ((thread_args_t*) args)->Edge2Adj;
	width = ((thread_args_t*) args)->width;
	height = ((thread_args_t*) args)->height;
	L = ((thread_args_t*) args)->L;
	l = ((thread_args_t*) args)->l;
	K = ((thread_args_t*) args)->K;
		
	for (k = lb; k <= ub; k++) 
		gridSteinerMT(g, Edge2Adj, width, height, T[k].terminals, T[k].no_terminals, L, T[k].net_num, &(T[k].edge_count), &(T[k].SteinerTree),l,K,H,A);
		
	return NULL;
}


#endif
