/*----------------------------------------------------------------------------------*/
/* graphtype.h																		*/
/*	Header file for grid graphs, and transform a grid graphs and virutal layer		*/
/*	grpahs.																			*/
/*																					*/
/*----------------------------------------------------------------------------------*/

#ifndef GRAPHTYPE_H_
#define GRAPHTYPE_H_

#include "graph.h"

int gen_gridgraph(Graph *grid_graph, int GRID_WIDTH, int GRID_HEIGHT);
void grid2VL(Graph *G, Graph *H);

void print_weights(Graph *g);
int find_edge_index( Graph *graph, CoordData *u, CoordData *v );
double find_edge_weight( Graph *graph, CoordData *u, CoordData *v );
int set_edge_weight( Graph *graph, CoordData *u, CoordData *v, int w);
void sync_weights(Graph *g, double *w);
void set_Edge2Adj(Graph *g, AdjList****a);
void set_capacity_lookup(int **c, Graph *g, int hc, int vc);
void set_length_lookup(double **c, Graph *g, double hc, double vc);

typedef struct graph_copies_args_t_ {
	Graph *grid;
	Graph *MLGraph;
} graph_copies_args_t;


void init_graph_copies_args(graph_copies_args_t **args, Graph *grid, Graph **MLGraph, int no_threads);
void *graph_runner(void *args);
void graph_copies(Graph *grid, Graph **H, int no_threads);



#endif
