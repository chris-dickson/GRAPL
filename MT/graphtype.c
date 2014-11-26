/*--------------------------------------------------------------------------------------*/
/* graphtype.c																			*/
/*	Implements routines to create grid graphs, and transform a grid graph to a			*/
/*  virutal layer grpah.																*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/


#ifndef GRAPHTYPE_C_
#define GRAPHTYPE_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "graph.h"
#include "list.h"
#include "set.h"
#include "graphalg.h"
#include "helper.h"

#include "graphtype.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 gen_gridgraph( ... )
		This function generates a grid graph with vertices [0..GRID_WIDTH) x [0..GRID_HEIGHT).
	Additionally, it adds 4 edges to each internal vertex. 2 edges to corner vertices, and 3 edges to
	each vertex on the boundary.   All edges added are given weight 1.   
	The graph is returned as the first parameter to the function.   There is no need to allocate/init
	a graph before calling this function
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int gen_gridgraph(Graph *grid_graph, int GRID_WIDTH, int GRID_HEIGHT) {

	PathVertex		*path_vertex, *path_start, *path_vertex1, *path_vertex2;
	CoordData		*coord;
	int				i,j;
	int				index;
	AdjList			***Ver2Adj2;		/* lookup to speed up edge insertion */
	
	printf("Creating Grid Graph\n");
	
	Ver2Adj2 = (AdjList***) malloc(sizeof(AdjList**)*(GRID_WIDTH));
	for (i = 0; i < GRID_WIDTH; i++)
		Ver2Adj2[i] = (AdjList**) malloc(sizeof(AdjList*)*GRID_HEIGHT);
	
	/* initialize a grid graph */
	graph_init(grid_graph, match_coord, (void*)path_vertex_free);
	printf("\tInitialized grid graph..\n");
					
	grid_graph->width = GRID_WIDTH;
	grid_graph->height = GRID_HEIGHT;
	grid_graph->layers = 1;				
	
    /* insert vertices into grid graph */
	printf("\tInserting verticies");
    for (i = 0; i < GRID_WIDTH; i++) {
		
		printf(".");
		fflush(stdout);
		
        for (j = 0; j < GRID_HEIGHT; j++) {
        
            /* allocate a vertex */
            if ((path_vertex = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
                printf("\t\tError allocating (%d,%d) PathVertex\nTerminating...\n",i,j);
                fflush(stdout);
                exit(1);
            }
                
            /* allocate a coordinate */
            if ((coord = (CoordData*) malloc(sizeof(CoordData))) == NULL) {
                printf("\t\tError allocating (%d,%d) CoordData\nTerminating...\n",i,j);
				fflush(stdout);
				exit(1);
            }
                
            /* set values to coord */
            coord->x = i;
            coord->y = j;
			coord->z = 0;
            
            /* insert coord as data to path vertex */
            path_vertex->data = coord;
			
			if ((i == 1)&&(j==2))
				path_start = path_vertex;
            
            /*insert it into the graph */
            if ((Ver2Adj2[i][j] = graph_ins_vertex2(grid_graph, path_vertex)) == NULL) {
                printf("\t\tError inserting vertex (%d,%d)\nTerminating...\n",i,j);
                return 1;
            }
        }
    } /*done inserting vertices*/     
	printf("Finished\n");
    
    
    /* for each vertex we check if there is a vertex to the left, right, above and below.
    In this way, it only consumes O(n) time where n is the number of vertices */
    /* insert edges into grid graph */
	printf("\tInserting edges");
	
	index = 1;
	
    for (i = 0; i < GRID_WIDTH; i++) {
		
		printf(".");
		fflush(stdout);
		
        for (j = 0; j < GRID_HEIGHT; j++) {
            
            /* check to the left */
            if (i - 1 >= 0) {
                CoordData *u,*v;
            
                /* allocate the vertices (for in the edges)*/
                u = (CoordData*)malloc(sizeof(CoordData));
                v = (CoordData*)malloc(sizeof(CoordData));
                
                if (u == NULL) {
                	printf("gsFH.h : u mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

                if (v == NULL) {
                	printf("gsFH.h : v mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

            
                /* set vertex data */
                u->x = i;
                u->y = j;
				u->z = 0;
				
                v->x = i-1;
                v->y = j;
				v->z = 0;
            
				
				if ((path_vertex1 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
                
                if ((path_vertex2 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
  					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
                    
                path_vertex1->data = u;
                path_vertex2->data = v;
                path_vertex2->weight = 1;
                
                if (graph_ins_edge2(grid_graph, path_vertex1, path_vertex2, Ver2Adj2) != 0)
                    return 1;
				
				free(path_vertex1->data);
				free(path_vertex1);
				
            }
            
            /* check to the right */
            if (i + 1 < GRID_WIDTH) {
                CoordData *u,*v;
            
                /* allocate the vertices (for in the edges)*/
                u = (CoordData*)malloc(sizeof(CoordData));
                v = (CoordData*)malloc(sizeof(CoordData));
                
                if (u == NULL) {
                	printf("gsFH.h : u mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

                if (v == NULL) {
                	printf("gsFH.h : v mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

            
                /* set vertex data */
                u->x = i;
                u->y = j;
				u->z = 0;
				
                v->x = i+1;
                v->y = j;
				v->z = 0;
            

				if ((path_vertex1 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
                
                if ((path_vertex2 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
  					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
                                    
                path_vertex1->data = u;
                path_vertex2->data = v;
                path_vertex2->weight = 1;
                
                if (graph_ins_edge2(grid_graph, path_vertex1, path_vertex2, Ver2Adj2) != 0)
                    return 1;
				
				free(path_vertex1->data);
				free(path_vertex1);
            }
    
            /* check up */
            if (j - 1 >= 0) {
                CoordData *u,*v;
            
                /* allocate the vertices (for in the edges)*/
                u = (CoordData*)malloc(sizeof(CoordData));
                v = (CoordData*)malloc(sizeof(CoordData));
                
                if (u == NULL) {
                	printf("gsFH.h : u mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

                if (v == NULL) {
                	printf("gsFH.h : v mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

            
                /* set vertex data */
                u->x = i;
                u->y = j;
				u->z = 0;
				
                v->x = i;
                v->y = j-1;
				v->z = 0;
            
			
				if ((path_vertex1 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
                
                if ((path_vertex2 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
  					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
                                    
                path_vertex1->data = u;
                path_vertex2->data = v;
                path_vertex2->weight = 1;
                
                if (graph_ins_edge2(grid_graph, path_vertex1, path_vertex2, Ver2Adj2) != 0)
                    return 1;

				free(path_vertex1->data);
				free(path_vertex1);

            }
            
            /* check down */
            if (j + 1 < GRID_HEIGHT) {
                CoordData *u,*v;
            
                /* allocate the vertices (for in the edges)*/
                u = (CoordData*)malloc(sizeof(CoordData));
                v = (CoordData*)malloc(sizeof(CoordData));

                if (u == NULL) {
                	printf("gsFH.h : u mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

                if (v == NULL) {
                	printf("gsFH.h : v mem allocation error\n");
                	fflush(stdout);
                	exit(1);
                }

            
                /* set vertex data */
                u->x = i;
                u->y = j;
				u->z = 0;
				
                v->x = i;
                v->y = j+1;
				v->z = 0;
            
				if ((path_vertex1 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
                
                if ((path_vertex2 = (PathVertex*) malloc(sizeof(PathVertex))) == NULL) {
  					printf("gsFH.h : path_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
                
                path_vertex1->data = u;
                path_vertex2->data = v;
                path_vertex2->weight = 1;
                
                if (graph_ins_edge2(grid_graph, path_vertex1, path_vertex2, Ver2Adj2) != 0)
                    return 1;
			
				free(path_vertex1->data);
				free(path_vertex1);

            }
    
    
        }
    }/* finished inserting edges */
	
	/* free lookup */
	for (i = 0; i < GRID_WIDTH; i++)
		free(Ver2Adj2[i]);
	free(Ver2Adj2);
	
	printf("Finished\n");
	return 0;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 grid2VL( ... )
		This function converts a grid graph G to a virtual (2 layer) graph H.  There is no need to 
	initialize/allocate H, but G should already be created.   If we are removing vertices (ie/ if the
	graph contains holes, this should be done BEFORE calling this function, as in this function, each 
	edge is enumerated with an integer.   Since we convert G to H, all edges in H are initially weight
	1.   Use sync_weights(...) to change the weights of edges after calling grid2VL().
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void grid2VL(Graph *G, Graph *H) { 
	ListElmt	*elmt;
	int			index,via,i,j,gw,gh,layers;
	AdjList		****Ver2Adj3;
	
	
	graph_init(H, match_coord, (void*)path_vertex_free);
	
	H->width = G->width;
	H->height = G->height;
	H->layers = 2;
	
	gw = H->width;
	gh = H->height;
	layers = H->layers;
	
	Ver2Adj3 = (AdjList****) malloc(sizeof(AdjList***)*gw);
	for (i = 0; i < gw; i++) {
		
		Ver2Adj3[i] = (AdjList***) malloc(sizeof(AdjList**)*gh);
		
		for (j = 0; j < gh; j++)
			Ver2Adj3[i][j] = (AdjList**) malloc(sizeof(AdjList*)*layers);
	}
	
	
	/* original edges enumrated from 1, vias enumerated after original edges.
	 also, since our graph is undirected, we have to half edges for each edge in undirected graph, 
	 thus we start vias from ecount/2 */
	index = 1;
	via = (G->ecount / 2) + 1;
	
	
	/* for each vertex in G, there will be two vertices in H.  (One per layer) */
	fflush(stdout);
	for (elmt = list_head(&G->adjlists); elmt != NULL; elmt = list_next(elmt)) {
		PathVertex	*v0,*v1,*v;		/* v = vertex from g, v1,v2 = new vertices to insert in H */
		CoordData	*v0c,*v1c;
				
		v = ((AdjList*)list_data(elmt))->vertex;
			
		/* insert both verticies */
		v0 = (PathVertex*) malloc(sizeof(PathVertex));
		v0c = (CoordData*) malloc(sizeof(CoordData));
		v1 = (PathVertex*) malloc(sizeof(PathVertex));
		v1c = (CoordData*) malloc(sizeof(CoordData));
		
		if ( (v0 == NULL)||(v0c == NULL)||(v1 == NULL)||v1c == NULL) {
			printf("gsFH.h : vertex mem allocation error (grid2VL)\n");
			fflush(stdout);
			exit(1);
		}
		
		v0->data = v0c;
		v1->data = v1c;
		
		/* horizontal vertex */
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;

		/* vertical vertex */
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		Ver2Adj3[v0c->x][v0c->y][v0c->z] = graph_ins_vertex2(H,v0);
		Ver2Adj3[v1c->x][v1c->y][v1c->z] = graph_ins_vertex2(H,v1);
		
		/* now we insert the via between v0 and v1 */
		v0 = (PathVertex*) malloc(sizeof(PathVertex));
		v0c = (CoordData*) malloc(sizeof(CoordData));
		v1 = (PathVertex*) malloc(sizeof(PathVertex));
		v1c = (CoordData*) malloc(sizeof(CoordData));

		if ( (v0 == NULL)||(v0c == NULL)||(v1 == NULL)||v1c == NULL) {
			printf("gsFH.h : vertex mem allocation error (grid2VL)\n");
			fflush(stdout);
			exit(1);
		}
		
		
		v0->data = v0c;
		v1->data = v1c;

		/* set the weights to 1, temporarily */
		v1->weight	 = 1;

		/* set via index */
		v1->index = via;

		/* horizontal vertex */
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;
		
		/* vertical vertex */
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		
		/* forward edge */
		graph_ins_edge3(H,v0,v1, Ver2Adj3);
		free(v0c);
		free(v0);


		v0 = (PathVertex*) malloc(sizeof(PathVertex));
		v0c = (CoordData*) malloc(sizeof(CoordData));
		v1 = (PathVertex*) malloc(sizeof(PathVertex));
		v1c = (CoordData*) malloc(sizeof(CoordData));
		
		if ( (v0 == NULL)||(v0c == NULL)||(v1 == NULL)||v1c == NULL) {
			printf("gsFH.h : vertex mem allocation error (grid2VL)\n");
			fflush(stdout);
			exit(1);
		}

		
		
		v0->data = v0c;
		v1->data = v1c;

		/* set the weights to 1, temporarily */
		v0->weight	 = 1;
		
		/* via index */
		v0->index = via;

		/* horizontal vertex */
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;
		
		/* vertical vertex */
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		/* forward edge */
		graph_ins_edge3(H,v1,v0, Ver2Adj3);
		free(v1c);
		free(v1);
		
		/* next via*/
		via++;
		
	}
	fflush(stdout);
	
	
	/* now we insert edges into the multi layer graph*/
	for (elmt = list_head(&G->adjlists); elmt != NULL; elmt = list_next(elmt)) {
		PathVertex	*v0,*v1,*v,*next;		/* v = vertex from g, v0,v1 = new vertices to insert in H*/
		CoordData	*v0c,*v1c;
		ListElmt	*e;
		int			vx,vy,vz;
		
		v = ((AdjList*)list_data(elmt))->vertex;
		
		vx = ((CoordData*)((PathVertex*)v)->data)->x;
		vy = ((CoordData*)((PathVertex*)v)->data)->y;
		vz = ((CoordData*)((PathVertex*)v)->data)->z;
		
		/* go through each vertex adjacent to v*/
		for (e = list_head(&((AdjList*)list_data(elmt))->adjacent); e != NULL; e = list_next(e)) {
			int			nx, ny, nz;
			double		weight;
			int			layer;
			int			exists;
			
			next = list_data(e);
			
			nx = ((CoordData*)((PathVertex*)next)->data)->x;
			ny = ((CoordData*)((PathVertex*)next)->data)->y;
			nz = ((CoordData*)((PathVertex*)next)->data)->z;
			weight = (double)((PathVertex*)next)->weight;

			/* find out if it's a horizontal or vertical edge*/
			if (vx == nx) layer = 1;
			else if (vy == ny) layer = 0;
			else assert(0);
			
			/* insert forward edge from v0 to v1*/
			v0 = (PathVertex*) malloc(sizeof(PathVertex));
			v1 = (PathVertex*) malloc(sizeof(PathVertex));
			v0c = (CoordData*) malloc(sizeof(CoordData));
			v1c = (CoordData*) malloc(sizeof(CoordData));
			
			if ( (v0 == NULL)||(v0c == NULL)||(v1 == NULL)||v1c == NULL) {
				printf("gsFH.h : vertex mem allocation error (grid2VL)\n");
				fflush(stdout);
				exit(1);
			}

			
			v0c->x = vx;
			v0c->y = vy;
			v0c->z = layer;
			
			v1c->x = nx;
			v1c->y = ny;
			v1c->z = layer;
			
			v0->data = v0c;
			v1->data = v1c;
			
			v1->weight = weight; 
			v1->index = index;
			
			graph_ins_edge3(H,v0,v1,Ver2Adj3);
			
			/* free the lookup variable*/
			free(v0c);
			free(v0);
			
			
			/* insert the backedge v1 to v0 */
			v0 = (PathVertex*) malloc(sizeof(PathVertex));
			v1 = (PathVertex*) malloc(sizeof(PathVertex));
			v0c = (CoordData*) malloc(sizeof(CoordData));
			v1c = (CoordData*) malloc(sizeof(CoordData));
			
			if ( (v0 == NULL)||(v0c == NULL)||(v1 == NULL)||v1c == NULL) {
				printf("gsFH.h : vertex mem allocation error (grid2VL)\n");
				fflush(stdout);
				exit(1);
			}

			
			v0c->x = vx;
			v0c->y = vy;
			v0c->z = layer;
			
			v1c->x = nx;
			v1c->y = ny;
			v1c->z = layer;
			
			v0->data = v0c;
			v1->data = v1c;
			
			v0->weight = weight;
			v0->index = index; 
			
			exists = graph_ins_edge3(H,v1,v0,Ver2Adj3);
			
			if (!exists) {
				((PathVertex*)next)->index = index;		/* set index in orginal graph G*/
				index++;
			}
			
			/* free the lookup variable*/
			free(v1c);
			free(v1);
			
		}
	}
	fflush(stdout);
	
	/* free lookup */
	for (i = 0; i < gw; i++) {
		
		for (j = 0; j < gh; j++)
			free(Ver2Adj3[i][j]);
			
		free(Ver2Adj3[i]);
	}
	free(Ver2Adj3);		
}


int find_edge_index( Graph *graph, CoordData *u, CoordData *v ) {
    ListElmt *e;
    AdjList *adjlistPtr;     
    
    adjlistPtr = NULL; 
    
    /* get adjlist of vertex associated with u */
    for (e = list_head(&graph->adjlists); e != NULL; e = list_next(e)) {
		int tx,ty,tz;
        AdjList *tmp;
		
        tmp = (AdjList*)list_data(e);
        tx = ((CoordData*)((PathVertex*)tmp->vertex)->data)->x;
        ty = ((CoordData*)((PathVertex*)tmp->vertex)->data)->y;
		tz = ((CoordData*)((PathVertex*)tmp->vertex)->data)->z;
        
        if ((u->x == tx)&&(u->y == ty)&&(u->z == tz))
            adjlistPtr = tmp;
    }
    
    if (adjlistPtr == NULL) {
        printf("Couldn't find vertex (%d,%d)'s adjacency list\n",u->x,u->y);
        return -1;
    }
    else {
        
        /* find vertex v in u's adjacency list */
        for (e = list_head(&adjlistPtr->adjacent); e != NULL; e = list_next(e)) {
            PathVertex *tmp;
            short int tx,ty,tz;
            
            tmp = list_data(e);
            
            tx = ((CoordData*)((PathVertex*)tmp)->data)->x;
            ty = ((CoordData*)((PathVertex*)tmp)->data)->y;
            tz = ((CoordData*)((PathVertex*)tmp)->data)->z;
            
            if ((tx == v->x)&&(ty == v->y)&&(tz == v->z)) 
                return ((PathVertex*)tmp)->index;
        }
		
        printf("Couldn't locate v in u's adjacency list\n");
        return -1;
    }
}




void *graph_runner(void *args)
{
	Graph *G, *H;
	
	G = ((graph_copies_args_t*) args)->grid;
	H = ((graph_copies_args_t*) args)->MLGraph;

	grid2VL(G,H);
	
	return NULL;
}

#endif
