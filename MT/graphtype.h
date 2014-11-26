/*--------------------------------------------------------------------------------------*/
/* graphtype.h																			*/
/*	Header file for grid graphs, and transform a grid graphs and virutal layer grpahs.	*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/

#ifndef GRAPHTYPE_H_
#define GRAPHTYPE_H_

#include "graph.h"
#include "helper.h"

int gen_gridgraph(Graph *grid_graph, int GRID_WIDTH, int GRID_HEIGHT);
void grid2VL(Graph *G, Graph *H);

int	find_edge_index( Graph *g, CoordData *u, CoordData *v );


typedef struct graph_copies_args_t_ {
	Graph *grid;
	Graph *MLGraph;
} graph_copies_args_t;

void *graph_runner(void *args);

#endif
