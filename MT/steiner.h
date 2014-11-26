/*--------------------------------------------------------------------------------------*/
/* steiner.h																			*/
/*	Header file of a 2-approximate steiner tree algorithm.   This is the				*/
/* block solver for our routing code.  The main idea is as follows:						*/
/*																						*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/


#ifndef STEINER_H_
#define STEINER_H_

#include "graph.h"
#include "helper.h"

/*--------------------------------------*/
/* unused in this code					*/
/*--------------------------------------*/
/* */ void create_spath_table(int m);/* */	
/* */ void reset_spath_table(int m); /* */
/* */ void destroy_spath_table();	 /* */
/*--------------------------------------*/

/* gridSteinerFH ---------------------------------------------------------------------------*/
/*	Returns : (int) 0 on success															*/
/*  Arguments :																				*/
/*		(Graph*)	  grid_graph: a pointer to the graph (not neccessarily the grid graph	*/
/*		(AdjList****) Edge2Adj	: 3d (x,y,z) array of pointers to adjacency lists in		*/
/*								grid_graph													*/
/*		(int) width,height		: width and height of grid_graph							*/
/*		(CoordData*) term		: array of terminals										*/
/*		(int) no_terminals		: number of terminals										*/
/*		(double*) L				: length function (not used)								*/
/*		(int) net_num			: index of the net in the instance							*/
/*		(int*) edge_count_out	: number of edges in the Steiner tree						*/
/*		(int**) SteinerTree_out	: pointer to an array of edges representing the Steiner tree*/
/*		(int) l					: current iteration of main algorithm						*/
/*		(int) K					: number of nets											*/
/*------------------------------------------------------------------------------------------*/
int gridSteinerFH(Graph *grid_graph, AdjList**** Edge2Adj, int width, int height, CoordData *term, 
	int no_terminals, double *L,int net_num, int *edge_count_OUT, int **SteinerTree_OUT, int l, int K);
	
#endif
