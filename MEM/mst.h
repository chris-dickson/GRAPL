/*--------------------------------------------------------------------------------------*/
/* mst.h																				*/
/*	Header file for an implementation of a prims algorithm for minimum spanning trees	*/
/*  This code was obtained from the book "Algorithms in C"								*/
/*																						*/
/*--------------------------------------------------------------------------------------*/



#ifndef MST_H_
#define MST_H_

#include "graph.h"
#include "graphalg.h"
#include "list.h"

int mst(Graph *graph, const MstVertex *start, List *span, int (*match)(const void *key1, const void *key2));


#endif
