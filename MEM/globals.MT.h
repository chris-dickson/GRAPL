/*----------------------------------------------------------------------------------*/
/*	globals.MT.h																	*/
/*		Header file for all global variables for MULTI-THREADED CODE				*/
/*----------------------------------------------------------------------------------*/

#ifndef GLOBALS_MT_H_
#define GLOBALS_MT_H_

#include <pthread.h>
#include "binheap.h"
#include "helper.h"
#include "routing.h"


/* Array of heaps and weight arrays.   Each thread gets its own heap and weight array.  */
Heap	**H;
double	**A_weights;

void	init_threaded_heaps(Heap ***H, double ***A, int gw, int gh, int no_threads);		
void	destroy_threaded_heaps(int no_threads);								



typedef struct SPElement{
	CoordData 	c;
	double	  	d;
	int		parent;
}SPElement;

int					**xy2term;							/* 2d array that maps terminals (x,y) position to an internal integer representation */
int					*computed_shortest_path;			/* uses internal int representation of a terminal to find if we have the shortest path estimates from that vertex */
SPElement			**shortest_paths;					
pthread_mutex_t		*computed_shortest_path_lock;			/* read/write locks for computed_shortest_path array (one lock for each array element) */
pthread_mutex_t		*shortest_paths_lock;

Net					*T;


/* debug */
int			useless_calls;						/* the number of wasted calls to dijkstra */


void	init_shortest_path_globals(int gw, int gh, int no_unique_terminals);
void	reset_shortest_path(int no_unique_terminals);
int		have_shortest_path(int x, int y);
void	found_shortest_path(int x, int y);		
void	init_xy2term(int gw, int gh, CoordData* unique_terminals, int no_unique_terminals);
void	destroy_shortest_path_globals(int gw, int gh, int no_unique_terminals);

typedef struct thread_args_t_{
	Heap *H;				/* Heap for the thread */
	double *A;				/* Weight array for the thread */
	int		lb,ub;			/* nets to be run in T, we are given lower and upper bounds inclusively */
	int		no_jobs;		/* number of nets to be run (ub - lb + 1) */
	
	Graph	*g;				/* pointer to the graph */
	AdjList ****Edge2Adj;	/* Adjacency list lookup */
	int		width;			/* width of g */
	int		height;			/* height of g */
	double	*L;				/* edge length lookup */
//	Net		*T;				/* net array */
	int		l;				/* iteration number */
	int		K;
}thread_args_t;

void	init_thread_args(thread_args_t **args, int no_threads, Heap **H, double **A, Graph **g, AdjList**** Edge2Adj, int width, int height, double *L,  int K, int l);
void	*runner(void *args);

#endif
