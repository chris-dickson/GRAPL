/*--------------------------------------------------------------------------------------*/
/* shortest.c																			*/
/*	An implementation of Dijkstra's algorithm for solving shortest path problems in		*/
/* graphs.   We use a pairing heap to speed up the process.   This is the bottleneck	*/
/* of the entire code.   If you make this faster, everything else will run faster too.	*/
/*  Portions of this code were obtained from the book "Algorithms in C"					*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/


#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "graph.h"
#include "graphalg.h"
#include "list.h"
#include "set.h"
#include "helper.h"
#include "binheap.h"



/*****************************************************************************
*                                                                            *
*  --------------------------------- relax --------------------------------  *
*                                                                            *
*****************************************************************************/

static int relax(PathVertex *u, PathVertex *v, double weight) {

/*****************************************************************************
*                                                                            *
*  Relax an edge between two vertices u and v.                               *
*                                                                            *
*****************************************************************************/

if (v->d > u->d + weight) {

   v->d = u->d + weight;
   v->parent = u;
   return 1;

}

return 0;

}

/*****************************************************************************
*                                                                            *
*  ------------------------------- shortest -------------------------------  *
*                                                                            *
*****************************************************************************/

int shortest(Graph *graph, const PathVertex *start, List *paths, int (*match)
   (const void *key1, const void *key2), int gw, int gh) {

AdjList            *adjlist;

PathVertex         *pth_vertex,
                   *adj_vertex;

ListElmt           *element,
                   *member;

int                found,
                   i;


/*PairHeap			H;
Position			*P;	*/

Heap				*H;					/* A_weights binary heap, used for priority queue */
double				*A_weights;			/* array of weights to be min-heapified */

CoordData			Index2Coord[(gw*gh*2)+1];
int					Coord2Index[gw][gh][2];

int					index;
AdjList				*Coord2Vertex[gw][gh][2];



/* initialize */
H = binheap_init(gw*gh*2);
A_weights = (double*) malloc(sizeof(double)*((gw*gh*2)+1));		
	


/*****************************************************************************
*                                                                            *
*  Initialize all of the vertices in the graph.                              *
*                                                                            *
*****************************************************************************/

found = 0;
index = 1;

for (element = list_head(&graph_adjlists(graph)); element != NULL; element =
   list_next(element)) {

   pth_vertex = ((AdjList *)list_data(element))->vertex;
   

   if (match(pth_vertex, start)) {
	   int x,y,z;

      /***********************************************************************
      *                                                                      *
      *  Initialize the start vertex.                                        *
      *                                                                      *
      ***********************************************************************/

	   x = ((CoordData*)((PathVertex*)pth_vertex)->data)->x;
	   y = ((CoordData*)((PathVertex*)pth_vertex)->data)->y;
	   z = ((CoordData*)((PathVertex*)pth_vertex)->data)->z;
	   Coord2Vertex[x][y][z] = list_data(element);
	   
	   
	  pth_vertex->color = white;
      pth_vertex->d = 0;
      pth_vertex->parent = NULL;
      found = 1;
	  

	  A_weights[index] = pth_vertex->d;

	  Coord2Index[x][y][z] = index;
	  Index2Coord[index].x = x;
	  Index2Coord[index].y = y;
	  Index2Coord[index].z = z;
	  index++;

      }

   else {
	   int x,y,z;

      /***********************************************************************
      *                                                                      *
      *  Initialize vertices other than the start vertex.                    *
      *                                                                      *
      ***********************************************************************/

	   x = ((CoordData*)((PathVertex*)pth_vertex)->data)->x;
	   y = ((CoordData*)((PathVertex*)pth_vertex)->data)->y;
	   z = ((CoordData*)((PathVertex*)pth_vertex)->data)->z;
	   Coord2Vertex[x][y][z] = list_data(element);

	   
	   
      pth_vertex->color = white;
      pth_vertex->d = DBL_MAX;
      pth_vertex->parent = NULL;

	  A_weights[index] = pth_vertex->d;
	  
	  Coord2Index[x][y][z] = index;
	  Index2Coord[index].x = x;
	  Index2Coord[index].y = y;
	  Index2Coord[index].z = z;
	  index++;


   }
}

/*****************************************************************************
*                                                                            *
*  Return if the start vertex was not found.                                 *
*                                                                            *
*****************************************************************************/

if (!found)
   return -1;


binheap_build(H,A_weights,gw*gh*2);		/* build the heap */

/*****************************************************************************
*                                                                            *
*  Use Dijkstra's algorithm to compute shortest paths from the start vertex. *
*                                                                            *
*****************************************************************************/

i = 0;

while (i < graph_vcount(graph)) {
	int x,y,z;

   /**************************************************************************
   *                                                                         *
   *  Select the white vertex with the smallest shortest-path estimate.      *
   *                                                                         *
   **************************************************************************/


	index = binheap_indexofmin(H);			/* get index of minimum-weight edge */
	binheap_extract(H);						/* remove it from the heap */
	
	x = Index2Coord[index].x;
	y = Index2Coord[index].y;
	z = Index2Coord[index].z;
	adjlist = Coord2Vertex[x][y][z];
	
	/**************************************************************************
   *                                                                         *
   *  Color the selected vertex black.                                       *
   *                                                                         *
   **************************************************************************/

   ((PathVertex *)adjlist->vertex)->color = black;
	

   /**************************************************************************
   *                                                                         *
   *  Traverse each vertex adjacent to the selected vertex.                  *
   *                                                                         *
   **************************************************************************/

   for (member = list_head(&adjlist->adjacent); member != NULL; member =
      list_next(member)) {
	   int		px,py,pz;

      adj_vertex = list_data(member);
	  
	  px = ((CoordData*)((PathVertex*)adj_vertex)->data)->x;
	  py = ((CoordData*)((PathVertex*)adj_vertex)->data)->y;
	  pz = ((CoordData*)((PathVertex*)adj_vertex)->data)->z;

      /***********************************************************************
      *                                                                      *
      *  Find the adjacent vertex in the list of adjacency-list structures.  *
      *                                                                      *
      ***********************************************************************/
	   
	  pth_vertex = ((AdjList*)Coord2Vertex[px][py][pz])->vertex;
					
	  /*****************************************************************
	  *                                                                *
	  *  Relax the adjacent vertex in the list of adjacency-list       *
	  *  structures.                                                   *
	  *                                                                *
	  *****************************************************************/

	  if (relax(adjlist->vertex, pth_vertex, adj_vertex->weight)) {
		  
		  /* update pth_vertex->d in FH */
		  binheap_decrease_key(H, Coord2Index[px][py][pz] , pth_vertex->d);
		}
   }

   /**************************************************************************
   *                                                                         *
   *  Prepare to select the next vertex.                                     *
   *                                                                         *
   **************************************************************************/

   i++;

}

/* destroy binary heap */
binheap_destroy(H);
free(A_weights);

/*****************************************************************************
*                                                                            *
*  Load the vertices with their path information into a list.                *
*                                                                            *
*****************************************************************************/


list_init(paths, NULL);

for (element = list_head(&graph_adjlists(graph)); element != NULL; element =
   list_next(element)) {

   /**************************************************************************
   *                                                                         *
   *  Load each black vertex from the list of adjacency-list structures.     *
   *                                                                         *
   **************************************************************************/

   pth_vertex = ((AdjList *)list_data(element))->vertex;

   if (pth_vertex->color == black) {

      if (list_ins_next(paths, list_tail(paths), pth_vertex) != 0) {
         printf("Problem inserting!\n");
         list_destroy(paths);
         return -1;
      }
   }
}

return 0;

}
