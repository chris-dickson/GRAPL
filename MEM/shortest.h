/*--------------------------------------------------------------------------------------*/
/* shortest.c																			*/
/*	Header for Dijkstra's algorithm for solving shortest path problems in				*/
/* graphs.																				*/
/*  Portions of this code were obtained from the book "Algorithms in C"					*/
/*	This file has been heavily modified since then...									*/
/*--------------------------------------------------------------------------------------*/


#ifndef SHORTEST_H_
#define SHORTEST_H_

#include "globals.h"
#include "graph.h"
#include "graphalg.h"
#include "list.h"

typedef struct sPathData_ {
	CoordData			*vertex,
	*parent;
	
	double				d;
} sPathData;

void *sPath_vertex_free(sPathData *s);
int copy_sPath(List *path, List *sPath, CoordData *start);			/* old way, create one shortest path at a time */
void create_path_list(List *dijkstra, SPElement *path);				/* new way, reusing output from dijkstra */
void SPElement2sPath(SPElement *path, List *sPath, CoordData *start );					/* convert the precomputed shortest path to an sPath */

int shortest(Graph *graph, const PathVertex *start, List *paths, int (*match) (const void *key1, const void *key2), int gw, int gh);


#endif
