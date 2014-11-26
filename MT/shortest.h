/*--------------------------------------------------------------------------------------*/
/* shortest.c																			*/
/*	Header for Dijkstra's algorithm for solving shortest path problems in				*/
/* graphs.   We use a pairing heap to speed up the process.   This is the bottleneck	*/
/* of the entire code.   If you make this faster, everything else will run faster too.	*/
/*  Portions of this code were obtained from the book "Algorithms in C"					*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/


#ifndef SHORTEST_H_
#define SHORTEST_H_

#include "graph.h"
#include "graphalg.h"
#include "list.h"




int shortest(Graph *graph, const PathVertex *start, List *paths, int (*match) (const void *key1, const void *key2), int gw, int gh);


#endif
