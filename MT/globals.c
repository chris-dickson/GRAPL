#ifndef GLOBALS_C_
#define GLOBALS_C_

#include <pthread.h>
#include <stdlib.h>
#include "globals.h"
#include "helper.h"
#include "fatal.h"
#include "steiner.h"
#include "graphtype.h"

void init_thread_args(thread_args_t **args, int no_threads, Graph **g, AdjList ****Edge2Adj, int width, int height, double *L, int K, int l)
{
	int i;
	
	(*args) = (thread_args_t*) malloc(sizeof(thread_args_t)*no_threads);
	
	for (i = 0; i < no_threads; i++)
	{
		thread_args_t *a;
		
		a = &((*args)[i]);
		
		a->g = g[i];
		a->Edge2Adj = Edge2Adj;
		a->width = width;
		a->height = height;
		a->L = L;
		a->l = l;
		a->K = K;
		
		a->lb = i*(K / no_threads) + 1;
		a->ub = (i+1)*(K / no_threads);
		
		if (i == no_threads - 1)
			a->ub += K % no_threads;
	}
}

void *runner(void *args) 
{
	int		lb,ub;
	Graph	*g;
	AdjList ****Edge2Adj;
	int		width,height;
	double	*L;
	int		l,K,k;
	
	lb = ((thread_args_t*) args)->lb;
	ub = ((thread_args_t*) args)->ub;
	g  = ((thread_args_t*) args)->g;
	Edge2Adj = ((thread_args_t*) args)->Edge2Adj;
	width = ((thread_args_t*) args)->width;
	height = ((thread_args_t*) args)->height;
	L = ((thread_args_t*) args)->L;
	l = ((thread_args_t*) args)->l;
	K = ((thread_args_t*) args)->K;
	
	for (k = lb; k <= ub; k++)
		if (T[k].fixed)
		{
			T[k].SteinerTree = T[k % K].SteinerTree;
			T[k].edge_count = T[k % K].edge_count;
		}
		else
			gridSteinerFH(g, Edge2Adj, width, height, T[k].terminals, T[k].no_terminals, L, T[k].net_num, &(T[k].edge_count), &(T[k].SteinerTree), l, K);
		
	return NULL;
}




#endif
