/* shortest.c																			*/
/*		Header file for computing shortest paths in graphs using Dijkstra's algorithm	*/
/*  Portions of this code were obtained from the book "Algorithms in C"					*/
/*																						*/
/*		This version is for MULTI-THREADED CODE											*/
/*--------------------------------------------------------------------------------------*/


#ifndef SHORTEST_MT_H_
#define SHORTEST_MT_H_

#include "globals.MT.h"
#include "graph.h"
#include "graphalg.h"
#include "list.h"

typedef struct sPathData_{
	CoordData			*vertex,*parent;
	double				d;
} sPathData;

void *sPath_vertex_free(sPathData *s);
int copy_sPath(List *path, List *sPath, CoordData *start);			/* old way, create one shortest path at a time */
void create_path_list(List *dijkstra, SPElement *path);				/* new way, reusing output from dijkstra */
void SPElement2sPath(SPElement *path, List *sPath, CoordData *start );					/* convert the precomputed shortest path to an sPath */

int shortest(Graph *graph, Heap *H, double *A_weights, const PathVertex *start, List *paths, int (*match) (const void *key1, const void *key2), int gw, int gh);


#endif
