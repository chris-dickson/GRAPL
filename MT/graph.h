/*--------------------------------------------------------------------------------------*/
/* graph.h																				*/
/*	Header for the basic graph routines													*/
/*  This code was obtained from the book "Algorithms in C"								*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>

#include "list.h"
#include "set.h"

/*****************************************************************************
*                                                                            *
*  Define a structure for adjacency lists.                                   *
*                                                                            *
*****************************************************************************/

typedef struct AdjList_ {

void               *vertex;
Set                adjacent;

} AdjList;

/*****************************************************************************
*                                                                            *
*  Define a structure for graphs.                                            *
*                                                                            *
*****************************************************************************/

typedef struct Graph_ {

int                vcount;
int                ecount;

int					width;
int					height;
int					layers;

int                (*match)(const void *key1, const void *key2);
void               (*destroy)(void *data);

List               adjlists;

} Graph;

/*****************************************************************************
*                                                                            *
*  Define colors for vertices in graphs.                                     *
*                                                                            *
*****************************************************************************/

typedef enum VertexColor_ {white, gray, black} VertexColor;

/*****************************************************************************
*                                                                            *
*  --------------------------- Public Interface ---------------------------  *
*                                                                            *
*****************************************************************************/

void graph_init(Graph *graph, int (*match)(const void *key1, const void
   *key2), void (*destroy)(void *data));

void graph_destroy(Graph *graph);

int graph_ins_vertex(Graph *graph, const void *data);
AdjList *graph_ins_vertex2(Graph *graph, const void *data);

int graph_ins_edge(Graph *graph, const void *data1, const void *data2);
int graph_ins_edge2(Graph *graph, const void *data1, const void *data2, AdjList ***Ver2Adj2);
int graph_ins_edge3(Graph *graph, const void *data1, const void *data2, AdjList ****Ver2Adj3);


int graph_rem_vertex(Graph *graph, void **data);

int graph_rem_edge(Graph *graph, void *data1, void **data2);

int graph_adjlist(const Graph *graph, const void *data, AdjList **adjlist);

int graph_is_adjacent(const Graph *graph, const void *data1, const void
   *data2);

#define graph_adjlists(graph) ((graph)->adjlists)

#define graph_vcount(graph) ((graph)->vcount)

#define graph_ecount(graph) ((graph)->ecount)

#endif
