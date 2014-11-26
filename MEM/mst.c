/*--------------------------------------------------------------------------------------*/
/* mst.h																				*/
/*	Implementation of prims algorithm for minimum spanning trees						*/
/*  This code was obtained from the book "Algorithms in C"								*/
/*																						*/
/*	This function computes the minimum spanning tree of the first argument graph.		*/
/*	It starts from the vertex given by start, and produces a list of MstVertices		*/
/*	into the list given by span.   The list is unordered and each vertex points			*/
/*	to its parent in the list.   The last argument is a function that is used			*/
/*	to match two mst_vertex structures.													*/
/*																						*/
/*--------------------------------------------------------------------------------------*/

#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include "graph.h"
#include "graphalg.h"
#include "list.h"

#include "mst.h"

int mst(Graph *graph, const MstVertex *start, List *span, int (*match)(const
   void *key1, const void *key2)) {

AdjList            *adjlist;

MstVertex          *mst_vertex,
                   *adj_vertex;

ListElmt           *element,
                   *member;

double             minimum;

int                found,
                   i;

/*****************************************************************************
*                                                                            *
*  Initialize all of the vertices in the graph.                              *
*                                                                            *
*****************************************************************************/

found = 0;
adjlist = NULL;

for (element = list_head(&graph_adjlists(graph)); element != NULL; element =
   list_next(element)) {

   mst_vertex = ((AdjList *)list_data(element))->vertex;

   if (match(mst_vertex, start)) {

      /***********************************************************************
      *                                                                      *
      *  Initialize the start vertex.                                        *
      *                                                                      *
      ***********************************************************************/

      mst_vertex->color = white;
      mst_vertex->key = 0;
      mst_vertex->parent = NULL;
      found = 1;

      }

   else {

      /***********************************************************************
      *                                                                      *
      *  Initialize vertices other than the start vertex.                    *
      *                                                                      *
      ***********************************************************************/

      mst_vertex->color = white;
      mst_vertex->key = DBL_MAX;
      mst_vertex->parent = NULL;

   }

}

/*****************************************************************************
*                                                                            *
*  Return if the start vertex was not found.                                 *
*                                                                            *
*****************************************************************************/

if (!found) {
    printf("Not found error\n");
   return -1;
}

/*****************************************************************************
*                                                                            *
*  Use Prim's algorithm to compute a minimum spanning tree.                  *
*                                                                            *
*****************************************************************************/

i = 0;

while (i < graph_vcount(graph)) {

   /**************************************************************************
   *                                                                         *
   *  Select the white vertex with the smallest key value.                   *
   *                                                                         *
   **************************************************************************/

   minimum = DBL_MAX;

   for (element = list_head(&graph_adjlists(graph)); element != NULL; element
      = list_next(element)) {

      mst_vertex = ((AdjList *)list_data(element))->vertex;

      if (mst_vertex->color == white && mst_vertex->key < minimum) {

         minimum = mst_vertex->key;
         adjlist = list_data(element);

      }

   }

   /**************************************************************************
   *                                                                         *
   *  Color the selected vertex black.                                       *
   *                                                                         *
   **************************************************************************/

   ((MstVertex *)adjlist->vertex)->color = black;

   /**************************************************************************
   *                                                                         *
   *  Traverse each vertex adjacent to the selected vertex.                  *
   *                                                                         *
   **************************************************************************/

   for (member = list_head(&adjlist->adjacent); member != NULL; member =
      list_next(member)) {

      adj_vertex = list_data(member);

      /***********************************************************************
      *                                                                      *
      *  Find the adjacent vertex in the list of adjacency-list structures.  *
      *                                                                      *
      ***********************************************************************/

      for (element = list_head(&graph_adjlists(graph)); element != NULL;
         element = list_next(element)) {

         mst_vertex = ((AdjList *)list_data(element))->vertex;

         if (match(mst_vertex, adj_vertex)) {

            /*****************************************************************
            *                                                                *
            *  Decide whether to change the key value and parent of the      *
            *  adjacent vertex in the list of adjacency-list structures.     *
            *                                                                *
            *****************************************************************/

            if (mst_vertex->color == white && adj_vertex->weight <
               mst_vertex->key) {

               mst_vertex->key = adj_vertex->weight;
               mst_vertex->parent = adjlist->vertex;

            }

            break;

         }

      }

   }

   /**************************************************************************
   *                                                                         *
   *  Prepare to select the next vertex.                                     *
   *                                                                         *
   **************************************************************************/

   i++;

}

/*****************************************************************************
*                                                                            *
*  Load the minimum spanning tree into a list.                               *
*                                                                            *
*****************************************************************************/

list_init(span, NULL);

for (element = list_head(&graph_adjlists(graph)); element != NULL; element =
   list_next(element)) {

   /**************************************************************************
   *                                                                         *
   *  Load each black vertex from the list of adjacency-list structures.     *
   *                                                                         *
   **************************************************************************/

   mst_vertex = ((AdjList *)list_data(element))->vertex;

   if (mst_vertex->color == black) {

      if (list_ins_next(span, list_tail(span), mst_vertex) != 0) {
         printf("Insertion error\n");
         list_destroy(span);
         return -1;

      }

   }

}

return 0;

}
