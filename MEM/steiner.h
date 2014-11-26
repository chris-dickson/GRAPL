/*--------------------------------------------------------------------------------------*/
/* steiner.h																			*/
/*	An implementation of Mehlhorn's 2-approximate steiner tree algorithm.   This is the	*/
/* block solver for our routing code (ie/   function evaluations are really steiner		*/
/* tree evaluations).  The main idea is as follows:										*/
/*--------------------------------------------------------------------------------------*/


#ifndef STEINER_H_
#define STEINER_H_

#include "graph.h"
#include "helper.h"

int gridSteinerFH(Graph *grid_graph, AdjList**** Edge2Adj, int width, int height, CoordData *term, 
	int no_terminals, double *L,int net_num, int *edge_count_OUT, int **SteinerTree_OUT, int l, int K);
	
#endif
