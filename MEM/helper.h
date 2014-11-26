/*--------------------------------------------------------------------------------------*/
/* helper.h																				*/
/*	Header for auxilliary routines that are used in many areas of code.					*/
/*																						*/
/*--------------------------------------------------------------------------------------*/


#ifndef HELPER_H_
#define HELPER_H_

#include "list.h"
#include "graph.h"
#include "graphalg.h"

void *path_vertex_free(PathVertex *p);
void *mst_vertex_free(MstVertex *m);

int match_coord( const void *grid1, const void *grid2);
int	match_coord2( const void *grid1, const void *grid2);

void my_itoa(int n, char *n_string);



#endif
