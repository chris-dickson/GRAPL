/*--------------------------------------------------------------------------------------*/
/* graphalg.h																			*/
/*	Defines structures for the graph algorithms (mst, dijkstra, etc)					*/
/*  This code was obtained from the book "Algorithms in C"								*/
/*--------------------------------------------------------------------------------------*/

#ifndef GRAPHALG_H
#define GRAPHALG_H

#include "graph.h"
#include "list.h"

/*****************************************************************************
*                                                                            *
*  Define a structure for vertices in minimum spanning trees.                *
*                                                                            *
*****************************************************************************/

typedef struct MstVertex_ {

void               *data;
double             weight;

VertexColor        color;
double             key;
int                is_leaf;
int				   index;		/* used to index edges */

struct MstVertex_  *parent;

} MstVertex;

/*****************************************************************************
*                                                                            *
*  Define a structure for vertices in shortest-path problems.                *
*                                                                            *
*****************************************************************************/

typedef struct PathVertex_ {

void               *data;
double             weight;

VertexColor        color;
double             d;
int				   is_leaf;
int				   index;		/* used to index edges */

struct PathVertex_ *parent;

} PathVertex;



/*****************************************************************************
*                                                                            *
*  --------------------------- Public Interface ---------------------------  *
*                                                                            *
*****************************************************************************/

int mst(Graph *graph, const MstVertex *start, List *span, int (*match)(const
   void *key1, const void *key2));

#endif
