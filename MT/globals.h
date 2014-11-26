#ifndef GLOBALS_H_
#define GLOBALS_H_

#include "graph.h"
#include "graphtype.h"
#include "routing.h"

Net *T;


typedef struct thread_args_t {

	int		lb,ub;
	Graph	*g;
	AdjList ****Edge2Adj;
	int		width;
	int		height;
	double	*L;
	int		l;
	int		K;
}thread_args_t;

void	init_thread_args(thread_args_t **args, int no_threads, Graph **g, AdjList ****Edge2Adj, int width, int height, double *L, int K, int l);
void	*runner(void *args);

#endif
