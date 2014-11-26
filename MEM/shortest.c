/* shortest.c																			*/
/*	An implementation of Dijkstra's algorithm for solving shortest path problems in		*/
/* graphs.   We use a binary heap to speed up the process.   This is the bottleneck		*/
/* of the entire code.   If you make this faster, everything else will run faster too.	*/
/*  Portions of this code were obtained from the book "Algorithms in C"					*/
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
#include "globals.h"
#include "shortest.h"


/*--------------------------------------------------------------*/
/*	Purpose: free sPath vertices								*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(sPathData *) s : vertex to be free						*/
/*--------------------------------------------------------------*/
void *sPath_vertex_free(sPathData *s) {
	free(s->vertex);
	free(s->parent);
	free(s);
	return(NULL);
}

/*--------------------------------------------------------------*/
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
/*						DEPRICATED								*/
/*--------------------------------------------------------------*/
/*	Purpose: make a copy of the shortest path returned by		*/
/*		Dijkstra's algorithm.									*/
/*  Returns: (List*) sPath										*/
/*  Args   :												    */
/*		(List*) path : Pointer to list returned by shortest()	*/
/*		(List*) sPath : pointer to list created					*/
/*		(CoordData*) start : start vertex in Dijkstra			*/
/*	Other Notes:												*/
/*			This is used so that the shortest paths don't get	*/
/*		overwritten.   The graph library stores the shortest	*/
/*		paths as pointers to the graph.   This function makes	*/
/*		a new list of sPath vertices that are independent of	*/
/*		of the graph.											*/
/*--------------------------------------------------------------*/
int copy_sPath(List *path, List *sPath, CoordData *start) {
	ListElmt		*element;
	sPathData		*newElement;
	PathVertex		*path_vertex, *v, *p;
	int				done;
	
	list_init(sPath, (void*)sPath_vertex_free);
	
	v = NULL;
	
	/* find start vertex in path */
	for ( element = list_head(path); element != NULL; element = list_next(element)) {
		int  sx,sy,sz,px,py,pz;
		
		path_vertex = list_data(element);
		px = ((CoordData*)((PathVertex*)path_vertex)->data)->x;
		py = ((CoordData*)((PathVertex*)path_vertex)->data)->y;	
		pz = ((CoordData*)((PathVertex*)path_vertex)->data)->z;	
		sx = start->x;
		sy = start->y;
		sz = start->z;
		
		
		/* if found, set a pointer and skip to end of list */
		if ((px==sx)&&(py==sy)&&(pz==sz)) {
			v = path_vertex;
			element = list_tail(path);
		}
	}
	
	p = v->parent;
	done = 0;	
	/* follow v back to root (when parent == NULL) */
	while ( !done ) {
		CoordData		*vCoord, *pCoord;
		
		/* allocate some memory */
		newElement = (sPathData*) malloc(sizeof(sPathData));
		vCoord = (CoordData*) malloc(sizeof(CoordData));
		pCoord = (CoordData*) malloc(sizeof(CoordData));
		
		if ((newElement == NULL)||(vCoord == NULL)||(pCoord == NULL)) {
			printf("helper.c : mem allocation error\n");
			fflush(stdout);
			exit(1);
		}
		
		/* set current vertex coords */
		vCoord->x = ((CoordData*)((PathVertex*)v)->data)->x;
		vCoord->y = ((CoordData*)((PathVertex*)v)->data)->y;
		vCoord->z = ((CoordData*)((PathVertex*)v)->data)->z;
		
		/* set parent vertex coords */
		pCoord->x = ((CoordData*)((PathVertex*)p)->data)->x;
		pCoord->y = ((CoordData*)((PathVertex*)p)->data)->y;
		pCoord->z = ((CoordData*)((PathVertex*)p)->data)->z;
		
		/* set the new element to list */
		newElement->vertex = vCoord;
		newElement->parent = pCoord;
		newElement->d = ((PathVertex*)v)->d;
		
		/* insert the element into the list */
		if (list_ins_next(sPath, list_tail(sPath), newElement) != 0) {
			printf("Problem inserting into sPath!\n");
			list_destroy(sPath);
			return -1;
		}
		
		
		/* next vertex */
		if (p->parent == NULL) {
			done = 1;
			
			/*insert one more element starting with parent */
			newElement = (sPathData*) malloc(sizeof(sPathData));
			vCoord = (CoordData*) malloc(sizeof(CoordData));
			
			if ((newElement==NULL)||(vCoord==NULL)) {
				printf("helper.c : memory allocation error\n");
				fflush(stdout);
				exit(1);
			}
			vCoord->x = ((CoordData*)((PathVertex*)p)->data)->x;
			vCoord->y = ((CoordData*)((PathVertex*)p)->data)->y;
			vCoord->z = ((CoordData*)((PathVertex*)p)->data)->z;
			newElement->vertex = vCoord;
			newElement->parent = NULL;
			newElement->d = 0;
			
			if (list_ins_next(sPath, list_tail(sPath), newElement) != 0){
				printf("Problem inserting into sPath!\n");
				list_destroy(sPath);
				return -1;
			}
		}
		else {
			v = p;
			p = p->parent;
		}
	}
	return 0;
}


/*--------------------------------------------------------------*/
/*	Purpose: Copy out shortest path data from the graph			*/
/*  Returns: (SPElement*) path									*/
/*  Args   :												    */
/*		(List*) dijkstra : list returned from dijkstra			*/
/*		(SPElement*) path : array of spelement vertices			*/
/*--------------------------------------------------------------*/
void create_path_list(List *dijkstra, SPElement *path) {
	ListElmt *e;
	PathVertex *p;
	
	for (e = list_head(dijkstra); e != NULL; e = list_next(e)) {
		
		p = (PathVertex*) list_data(e);
		
		path[ p->index ].c.x = ((CoordData*)(p->data))->x;
		path[ p->index ].c.y = ((CoordData*)(p->data))->y;
		path[ p->index ].c.z = ((CoordData*)(p->data))->z;
		path[ p->index ].d   = p->d;
		
		if (p->parent != NULL)
			path[ p->index ].parent = ((PathVertex*)(p->parent))->index;
		else
			path[ p->index ].parent = -1;
	}
}

/*--------------------------------------------------------------*/
/*	Purpose: Conrvet shortest path list to an sPath list 		*/
/*  Returns: (List*) sPath										*/
/*  Args   :												    */
/*		(SPElement*) path : path created by create_path_list()	*/
/*		(List*) sPath : list of sPath Vertices					*/
/*		(CoordData*) start : start vertex in list				*/
/*--------------------------------------------------------------*/
void SPElement2sPath(SPElement *path, List *sPath, CoordData *start) {
	int			s,p,v;			/* "pointers" to vertices in path */
	sPathData 	*newElement;	
	CoordData	*vCoord,*pCoord;
	int 		done;
	
	
	list_init(sPath, (void*)sPath_vertex_free);
	
	/* find end vertex in path (root of list, parent == -1) */
	s = 0;
	while ( path[s].parent != -1 )
		s++;
	
	
	/* find startd vertex in path */
	v = 0;
	while (!(( path[v].c.x == start->x )&&( path[v].c.y == start->y )&&( path[v].c.z == start->z )))
		v++;
		
	p = path[v].parent;
	
	done = 0;	
	/* follow v back to root (when v == s) */
	while (!done) {	
		
		/* allocate some memory */
		newElement = (sPathData*) malloc(sizeof(sPathData));
		vCoord = (CoordData*) malloc(sizeof(CoordData));
		pCoord = (CoordData*) malloc(sizeof(CoordData));
		
		if ((newElement == NULL)||(vCoord == NULL)||(pCoord == NULL)) {
			printf("shortest.c : mem allocation error\n");
			fflush(stdout);
			exit(1);
		}
		
		/* set current vertex coords */
		vCoord->x = path[v].c.x;
		vCoord->y = path[v].c.y;
		vCoord->z = path[v].c.z;
		
		/* set parent vertex coords */
		pCoord->x = path[p].c.x;
		pCoord->y = path[p].c.y;
		pCoord->z = path[p].c.z;
		
		/* set the new element to list */
		newElement->vertex = vCoord;
		newElement->parent = pCoord;
		newElement->d = path[v].d;
		
		/* insert the element into the list */
		if (list_ins_next(sPath, list_tail(sPath), newElement) != 0) {
			printf("Problem inserting into sPath!\n");
			list_destroy(sPath);
			exit(1);
		}
		
		
		/* next vertex */
		if (path[p].parent == -1) {
			done = 1;
			
			/*insert one more element starting with parent */
			newElement = (sPathData*) malloc(sizeof(sPathData));
			vCoord = (CoordData*) malloc(sizeof(CoordData));
			
			if ((newElement==NULL)||(vCoord==NULL)) {
				printf("helper.c : memory allocation error\n");
				fflush(stdout);
				exit(1);
			}
			vCoord->x = path[p].c.x;
			vCoord->y = path[p].c.y;
			vCoord->z = path[p].c.z;
			newElement->vertex = vCoord;
			newElement->parent = NULL;
			newElement->d = 0;
			
			if (list_ins_next(sPath, list_tail(sPath), newElement) != 0){
				printf("Problem inserting into sPath!\n");
				list_destroy(sPath);
				exit(1);
			}
		}
		else {
			v = p;
			p = path[p].parent;
		}
	}
	
}
	



/*--------------------------------------------------------------*/
/*	Purpose: relaxation routine called by dijkstra				*/
/*  Returns: 1 if vertices need to be relaxed, 0 otherwise		*/
/*  Args   :												    */
/*		(PathVertex*) u,v : source and destination				*/
/*		(double) weight : weight added by dijkstra				*/
/*--------------------------------------------------------------*/
static int relax(PathVertex *u, PathVertex *v, double weight) {

if (v->d > u->d + weight) {

   v->d = u->d + weight;
   v->parent = u;
   return 1;

}

return 0;

}

/*--------------------------------------------------------------*/
/*	Purpose: Dijkstra's algorithm for computing shortest paths	*/
/*  Returns: 0 on success.  (List*) paths						*/
/*  Args   :												    */
/*		(Graph*) graph : pointer to the graph					*/
/*		(const PathVertex*) start : start/source vertex			*/
/*		(List*) paths : List of PathVertex* to graph vertices	*/
/*		(int*) match : function pointer used to match vertices	*/
/*		(int) gw,gh : width and height of graph					*/
/*--------------------------------------------------------------*/
int shortest(Graph *graph, const PathVertex *start, List *paths, int (*match)
   (const void *key1, const void *key2), int gw, int gh) {

AdjList            *adjlist;

PathVertex         *pth_vertex,
                   *adj_vertex;

ListElmt           *element,
                   *member;

int                found,
                   i;

CoordData			Index2Coord[(gw*gh*2)+1];
int					Coord2Index[gw][gh][2];

int					index;
AdjList				*Coord2Vertex[gw][gh][2];


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
	   short   int x,y,z;

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
	   short   int x,y,z;

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
	short   int x,y,z;

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
	   short   int		px,py,pz;

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
/*binheap_destroy(H);*/
/*free(A_weights);*/

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
