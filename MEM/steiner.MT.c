/*--------------------------------------------------------------------------------------*/
/* steiner.c																			*/
/*	An implementation of Mehlhorn's 2-approximate steiner tree algorithm.   This is the	*/
/* block solver for our routing code (ie/   function evaluations are really steiner		*/
/* tree evaluations).																	*/
/*--------------------------------------------------------------------------------------*/

#ifndef STEINER_MT_C_
#define STEINER_MT_C_

#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

#include "graph.h"
#include "graphalg.h"
#include "graphtype.h"
#include "mst.h"
#include "shortest.MT.h"
#include "list.h"
#include "set.h"
#include "routing.h"

#include "globals.MT.h"
#include "steiner.MT.h"

#define PRINT printf("Here\n");fflush(stdout);
 
/*--------------------------------------------------------------*/
/*	Purpose: Function for computing one Steiner tree			*/
/*  Returns: 0 on success.  (int*) *SteinerTree_OUT				*/
/*  Args   :												    */
/*		(Graph*) grid_graph : pointer to the graph				*/
/*		(AdjList****) : adjacency list lookup function			*/
/*		(int) width,height : width and height of grid_graph		*/
/*		(CoordData*) term : array of terminals					*/
/*		(int) no_terminals : length of term array				*/
/*		(double*) L : edge length lookup						*/
/*		(int) net_num : identifier								*/
/*		(int*) edge_count_OUT : number of edges in Steiner tree	*/
/*					including vias								*/
/*		(int*) *Steiner_Tree_OUT : pointer to out Steiner tree	*/
/*					given as an array of integers representing	*/
/*					edge indices								*/
/*		(int) l : iteration number								*/
/*		(int) K : number of nets								*/
/*		(Heap*) H : pointer to the global heap to be used		*/
/*		(double*) A_weights : pointer to the weight array		*/
/*	Other Notes:												*/
/*		We provide an outline of the approximation algorithm.	*/
/*		The main idea is as follows:							*/
/*		1.   Generate complete distance network (ie/ shortest	*/
/*				path from each terminals to every other terminal*/ 
/*				in the graph.  (BOTTLENECK!!!)					*/
/*		2.   Compute minimum spanning tree of result from 1.	*/
/*		3.   Replace each edge of mst from 2 with the original	*/
/*				grid edges from input graph.					*/
/*		4.   Compute minimum spanning tree of result from 3.	*/
/*		5.   Remove each non-terminal leaf from the spanning	*/
/*				tree in 4.										*/
/*		We are left with a steiner tree that connects the		*/
/*		terminals specified.									*/
/*--------------------------------------------------------------*/
int gridSteinerMT(Graph *grid_graph, AdjList**** Edge2Adj, int width, int height, CoordData *term, int no_terminals, double *L,int
	net_num, int *edge_count_OUT, int **SteinerTree_OUT, int l, int K, Heap *H, double *A_weights) {
    int i,j;
	/* globals */
	
	int				GRID_WIDTH, GRID_HEIGHT, NO_TERMINALS;
	
	
	Graph			ND;					/* distance network of terminals */
	
	/*used to generate grid graph */
	PathVertex      *path_vertex;
	
	sPathData		*sPath_vertex;
	
	CoordData       *coord;
	
	/* used in shortest paths computations */
	PathVertex      **terminals;        /* a (1D) array of terminal (pointers) */
	List            *tempPaths,			/* tempporary */
					**sPaths;			/* 2d array of shortest path lists */
	ListElmt        *element;           /* temp, used to traverse a list */	
	
	/* used in timing */
	
	/* used in computing MST of ND */
	MstVertex		*mst_start,
					*mst_vertex,
					*mst_vertex1,
					*mst_vertex2;
	List			TD;					/* spanning tree of ND */
										/* used in final steps of algorithm */
	Graph			NTD;
	List			T;					/* spanning tree of NTD (eventually the steiner tree)*/
	int				isSteiner;			/* boolean flag */
	


/*------------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------------------------*/
	
			
	GRID_WIDTH = width;
	GRID_HEIGHT = height;
	NO_TERMINALS = no_terminals;
	
	mst_start = NULL;



	/* Create Path Vertices out of the original terminal set */
	if ((terminals = (PathVertex**) malloc(sizeof(PathVertex*)*NO_TERMINALS))==NULL) {
		printf("gsFH.h : terminals mem allocation error\n");
		fflush(stdout);
		exit(1);
	}	
	for (i = 0; i < NO_TERMINALS; i++) {
		short int x,y,z;
		
		x = term[i].x;
		y = term[i].y;
		z = term[i].z;
		
		path_vertex = (PathVertex*) malloc(sizeof(PathVertex));
		coord = (CoordData*) malloc(sizeof(CoordData));
		
		if ((path_vertex == NULL)||(coord == NULL)) {
			printf("gsFH.h : terminal[i] mem allocation error\n");
			fflush(stdout);
			exit(1);
		}

		
		coord->x = x;
		coord->y = y;
		coord->z = z;
		path_vertex->data = coord;
		
		terminals[i] = path_vertex;
	}

    

	/* inialize an array of list pointers used in extracting shortest paths from Dijkstra */
	sPaths = (List**) malloc(sizeof(List*)*NO_TERMINALS);
	if (sPaths == NULL) {
		printf("gsFH.h : sPaths mem allocation error\n");
		fflush(stdout);
		exit(1);
	}
	if ((tempPaths = (List*) malloc(sizeof(List)*NO_TERMINALS))==NULL) {
		printf("gsFH.h : tempPaths mem allocation error\n");
		fflush(stdout);
		exit(1);
	}	
	for (i = 0; i < NO_TERMINALS; i++)
		if ((sPaths[i] = (List*) malloc(sizeof(List)*NO_TERMINALS))==NULL) {
			printf("gsFH.h : sPaths[i] mem allocation error\n");
			fflush(stdout);
			exit(1);
		}



    
    /*--------------------------------------------------------------------------------------*/
    /*                          COMPUTE THE SHORTEST PATHS                                  */
    /*--------------------------------------------------------------------------------------*/
    
	
    
    /* for each terminal (stored as a path vertex), find the shortest path   */
    /* The shortest path for terminal[i] is stored in the List pointed to by */
    /* paths[i].                                                             */

	for (i = 0; i < NO_TERMINALS; i++) {
		
		/* lock this terminal */
		pthread_mutex_lock(&computed_shortest_path_lock[xy2term[term[i].x][term[i].y]]);		
		if (!have_shortest_path(term[i].x,term[i].y)) {
		
			if (shortest(grid_graph, H, A_weights, terminals[i], &tempPaths[i], match_coord,GRID_WIDTH,GRID_HEIGHT) != 0)
				return 1;		
		
			pthread_mutex_lock(&shortest_paths_lock[xy2term[term[i].x][term[i].y]]);
			create_path_list(&tempPaths[i], shortest_paths[ xy2term[ term[i].x ][ term[i].y ]]);
			pthread_mutex_unlock(&shortest_paths_lock[xy2term[term[i].x][term[i].y]]);
			found_shortest_path( term[i].x, term[i].y );
//			list_destroy(&tempPaths[i]);			
		}
		/* unlock it */
		pthread_mutex_unlock(&computed_shortest_path_lock[ xy2term[term[i].x][term[i].y] ]);
		
		
		/* copy out the shortest path data, if we don't do this, it will get overwritten*/
		for (j = 0; j < NO_TERMINALS; j++) 
			if (i != j) {
				pthread_mutex_lock(&shortest_paths_lock[xy2term[term[i].x][term[i].y]]);
				SPElement2sPath( shortest_paths[ xy2term[term[i].x][term[i].y] ], &(sPaths[i][j]), (CoordData*)((PathVertex*)terminals[j])->data); 
				pthread_mutex_unlock(&shortest_paths_lock[xy2term[term[i].x][term[i].y]]);
			}
			
	}

	free(tempPaths);
		
		
	/*--------------------------------------------------------------------------------------------------*/
	/*                Generate complete distance network ND												*/
	/*--------------------------------------------------------------------------------------------------*/

	/* initialize the graph */
	graph_init(&ND, match_coord, (void*)mst_vertex_free);
	
	/* insert the verticies.  Verticies consist of all the terminals */
	for (i = 0; i < NO_TERMINALS; i++) {
		
		/* allocate space for a MST vertex */
		if ((mst_vertex = (MstVertex*) malloc(sizeof(MstVertex))) == NULL) {
			printf("Error allocating space for mst_vertex\n");
			printf("Terminating..\n");
			return 1;
		}
		
		/* if it's the first, make it the start.  It doesn't matter which one is the start */
		if (i == 1)
			mst_start = mst_vertex;
		
		/* set the data */
		if ((coord = (CoordData*) malloc(sizeof(CoordData)))==NULL) {
			printf("gsFH.h : coord mem allocation error\n");
			fflush(stdout);
			exit(1);
		}

		coord->x = ((CoordData*)(((PathVertex*)terminals[i])->data))->x;
		coord->y = ((CoordData*)(((PathVertex*)terminals[i])->data))->y;
		coord->z = ((CoordData*)(((PathVertex*)terminals[i])->data))->z;
		mst_vertex->data = coord;
		
		/* insert */
		if (graph_ins_vertex(&ND, mst_vertex) != 0) {
			printf("Error inserting vertex into mst_graph\n");
			printf("Terminating...\n");
			return 1;
		}
	}

	/*	now we must insert the edges into the distance network graph ND.  We do this by accessing the
		shortest path lists (sPath) computed in the previous step */
	for (i = 0; i < NO_TERMINALS; i++) {
		short int ux,uy,uz;
		
		ux = ((CoordData*)((PathVertex*)terminals[i])->data)->x;
		uy = ((CoordData*)((PathVertex*)terminals[i])->data)->y;		
		uz = ((CoordData*)((PathVertex*)terminals[i])->data)->z;
				
		for (j = 0; j < NO_TERMINALS; j++) {
			short int vx,vy,vz;
			
			/* shouldn't be an edge from a terminal to itself */
			if (i != j) {
				double weight;
				CoordData *v1,*v2;
				int eCode;
			
				vx = ((CoordData*)((PathVertex*)terminals[j])->data)->x;
				vy = ((CoordData*)((PathVertex*)terminals[j])->data)->y;
				vz = ((CoordData*)((PathVertex*)terminals[j])->data)->z;
				
			
				/* now we must find how far away vx is from vy.   we do this by looking
					for at the head element in sPath[i][j] */
				
				element = list_head(&(sPaths[i][j]));
				sPath_vertex = list_data(element);
				weight = ((sPathData*)sPath_vertex)->d;
				
				/* allocate an edge */				
				if ((v1 = (CoordData*) malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : v1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				if ((v2 = (CoordData*) malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : v2 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				
				v1->x = ux;
				v1->y = uy;
				v1->z = uz;
				
				v2->x = vx;
				v2->y = vy;
				v2->z = vz;
				
				if ((mst_vertex1 = (MstVertex*) malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
								
				if ((mst_vertex2 = (MstVertex*) malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_vertex2 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				mst_vertex1->data =v1;
				mst_vertex2->data = v2;
				mst_vertex2->weight = weight;
				
				if ((eCode = graph_ins_edge(&ND, mst_vertex1, mst_vertex2)) != 0) {
					printf("Error inserting edge into ND\n");
					printf("graph_ins_edge returned the value %d\n",eCode);
					return 1;
				}
				
				free(mst_vertex1->data);
				free(mst_vertex1);
				
			}/* endif i!=j */
		}/* endfor j */
	}/* endfor i */

/*--------------------------------------------------------------------------------------------------*/
/*		Copmute TD (Min Span Tree of ND)															*/
/*--------------------------------------------------------------------------------------------------*/

	if (mst(&ND, mst_start,&TD, match_coord) != 0) {
		printf("Error computing minimum spanning tree\n");
		return 1;
	}	
	
	/* set leaves */
	
	/* initialize */
	for ( element = list_head(&TD); element != NULL; element = list_next(element))
	{
	   mst_vertex = list_data(element);
	   mst_vertex->is_leaf = 1;
    }
    
    /* for each node, set the parent is_leaf to 0.   Then, all leaves will remain */
    for (element = list_head(&TD); element != NULL; element = list_next(element))
    {
        mst_vertex = list_data(element);
        
        if (mst_vertex->parent != NULL)
            mst_vertex->parent->is_leaf = 0;
    }
	

/*--------------------------------------------------------------------------------------------------*/
/*	         Find N[TD]																				*/
/*--------------------------------------------------------------------------------------------------*/
	
	graph_init(&NTD,match_coord,(void*)mst_vertex_free);


    /* for each edge in TD */
	for (element = list_head(&TD); element != NULL; element = list_next(element)) {
		MstVertex	*nextVertex;
		int			v,p;
		

		nextVertex = list_data(element);
		
		/* if it is not the root */
		if (nextVertex->parent != NULL) {
			short   int vx,vy,vz,px,py,pz;
			ListElmt	*currentV, *nextV;
			int			done;
			
			
			vx = ((CoordData*)((MstVertex*)nextVertex)->data)->x;
			vy = ((CoordData*)((MstVertex*)nextVertex)->data)->y;
			vz = ((CoordData*)((MstVertex*)nextVertex)->data)->z;
			
			px = ((CoordData*)(MstVertex*)(nextVertex->parent)->data)->x;
			py = ((CoordData*)(MstVertex*)(nextVertex->parent)->data)->y;			
			pz = ((CoordData*)(MstVertex*)(nextVertex->parent)->data)->z;
			
			/* find terminal index of nextVertex and nextVertex->parent */
			for (i = 0; i < NO_TERMINALS; i++) {
				int tx,ty,tz;
				
				tx = ((CoordData*)((PathVertex*)terminals[i])->data)->x;
				ty = ((CoordData*)((PathVertex*)terminals[i])->data)->y;
				tz = ((CoordData*)((PathVertex*)terminals[i])->data)->z;
				
				if ((tx == vx)&&(ty == vy)&&(tz == vz))
					v = i;
				
				if ((tx == px)&&(ty == py)&&(tz == pz))
					p = i;
			}
			
			/* now, we must step through the list of sPathData elements found in sPaths[p][v].
			   For each element in the list, we must insert vertices for the vertex and parent, 
				then make an edge with the appropriate weight and insert it */
			
			currentV = list_head(&(sPaths[p][v]));
			nextV = list_next(currentV);
			done = 0;
			while ( !done ) {
				MstVertex	*u,*v, *mst_vertex1, *mst_vertex2;
				CoordData	*uc,*vc;
				sPathData	*currentVData, *nextVData;
				short   int			cvx,cvy,cvz,nvx,nvy,nvz;
				double		weight;

				/*---------------------------------------*/
				/* insert vertices u and v into NTD      */
				/*---------------------------------------*/
				
				/* make a vertex for currentV and nextV */
				u = (MstVertex*) malloc(sizeof(MstVertex));
				v = (MstVertex*) malloc(sizeof(MstVertex));
				uc = (CoordData*) malloc(sizeof(CoordData));
				vc = (CoordData*) malloc(sizeof(CoordData));
				
				if ((u == NULL)||(uc==NULL)||(v==NULL)||(vc==NULL)) {
					printf("gsFH.h : error allocating vertex for NTD\n");
					fflush(stdout);
					exit(1);
				}
				
				/* get vertices from the sPaths list */
				currentVData = list_data(currentV);
				nextVData = list_data(nextV);
				
				/* get vertex data */
				cvx = ((CoordData*)((sPathData*)currentVData)->vertex)->x;
				cvy = ((CoordData*)((sPathData*)currentVData)->vertex)->y;
				cvz = ((CoordData*)((sPathData*)currentVData)->vertex)->z;
				
				nvx = ((CoordData*)((sPathData*)nextVData)->vertex)->x;
				nvy = ((CoordData*)((sPathData*)nextVData)->vertex)->y;
				nvz = ((CoordData*)((sPathData*)nextVData)->vertex)->z;

								
				/* set vertex data */
				uc->x = cvx;
				uc->y = cvy;
				uc->z = cvz;
				
				vc->x = nvx;
				vc->y = nvy;
				vc->z = nvz;
				
				u->data = uc;
				v->data = vc;
				
				/* calculate weight between u and v */
				weight = currentVData->d - nextVData->d;

				
				/* try and insert u, if it exists, then delete the memory we allocated for it */
				if ( graph_ins_vertex(&NTD, u) == 1 ) {
					free(uc);
					free(u);
				}
				else {
					/* doesnt' matter which one is the start */
					mst_start = u;
				}

				/* try and insert v, if it exists, then delete the memorr we allocated for it */
				if ( graph_ins_vertex(&NTD, v) == 1) {
					free(vc);
					free(v);
				}
				
				/* now the vertices u and v are in the graph.   we now have to make an edge for uv */
				
				/* make edge going forward */				
				if ((uc = (CoordData*)malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : uc mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				if ((vc = (CoordData*)malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : vc mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				uc->x = cvx;
				uc->y = cvy;
				uc->z = cvz;

				vc->x = nvx;
				vc->y = nvy;
				vc->z = nvz;
				
				if ((mst_vertex1 = (MstVertex*) malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_Vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				if ((mst_vertex2 = (MstVertex*)malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_vertex2 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				mst_vertex1->data = uc;
				mst_vertex2->data = vc;
				mst_vertex2->weight = weight;
				
				/* try and insert, if it exists, free previously allocated mem */
				if ( graph_ins_edge(&NTD, mst_vertex1, mst_vertex2) == 1) {
					free(vc);
					free(mst_vertex2);
				}
				
				/* free the label */
				free(mst_vertex1->data);
				free(mst_vertex1);
				
				
				/* make edge going backward */				
				if ((uc = (CoordData*)malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : uc mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				if ((vc = (CoordData*)malloc(sizeof(CoordData))) == NULL) {
					printf("gsFH.h : vc mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				uc->x = cvx;
				uc->y = cvy;
				uc->z = cvz;
				
				vc->x = nvx;
				vc->y = nvy;
				vc->z = nvz;
				
				if ((mst_vertex1 = (MstVertex*) malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_Vertex1 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}

				if ((mst_vertex2 = (MstVertex*)malloc(sizeof(MstVertex))) == NULL) {
					printf("gsFH.h : mst_Vertex2 mem allocation error\n");
					fflush(stdout);
					exit(1);
				}
				
				mst_vertex1->data = vc;
				mst_vertex2->data = uc;
				mst_vertex2->weight = weight;
				
				/* try and insert, if it exists, free previously allocated mem */
				if ( graph_ins_edge(&NTD, mst_vertex1, mst_vertex2) == 1) {
					free(uc);
					free(mst_vertex2);
				}

				/* free the label */
				free(mst_vertex1->data);
				free(mst_vertex1);
								
				
				/* follow pointers */
				currentV = list_next(currentV);
				nextV = list_next(nextV);
				
				/* check to see if we're finished */
				if (nextV == NULL)
					done = 1;
			
			}
				
		 }
	}


/*----------------------------------------------------------------------------------------------------------*/
/*	Compute T (minimum spanning tree of NTD)																*/
/*----------------------------------------------------------------------------------------------------------*/

	/* call minimum spanning tree subroutine */
	if (mst(&NTD, mst_start,&T, match_coord) != 0) {
		printf("Error computing minimum spanning tree\n");
		return 1;
	}

	/* set leaves */	
	for ( element = list_head(&T); element != NULL; element = list_next(element))
	{
		mst_vertex = list_data(element);
		mst_vertex->is_leaf = 1;
    }
    
    /* for each node, set the parent is_leaf to 0.   Then, all leaves will remain */
    for (element = list_head(&T); element != NULL; element = list_next(element))
    {
        mst_vertex = list_data(element);
        
        if (mst_vertex->parent != NULL)
            mst_vertex->parent->is_leaf = 0;
    }
    
	
	
/*--------------------------------------------------------------------------------------*/
/*		Compute Steiner Tree Sk															*/
/*--------------------------------------------------------------------------------------*/

	isSteiner = 0;
	/* we remove all leaves that arent' terminals, when all leaves are terminals then we have a Steiner Tree */
	while (!isSteiner) {
		ListElmt	*prev;
		
		/* assume we have it */
		isSteiner = 1;
		
		/* check if each leaf is a terminal */
		prev = list_head(&T);
		element = list_next(prev);
		while (element != NULL) {
			short   int mx,my,mz;
			
			mst_vertex = list_data(element);
			mx = ((CoordData*)((MstVertex*)mst_vertex)->data)->x;
			my = ((CoordData*)((MstVertex*)mst_vertex)->data)->y;		
			mz = ((CoordData*)((MstVertex*)mst_vertex)->data)->z;					
			
			if (mst_vertex->is_leaf) {
				int found;
				
				found = 0;
				for (i = 0; i < NO_TERMINALS; i++) {
					int	tx,ty,tz;
					
					tx = ((CoordData*)((PathVertex*)terminals[i])->data)->x;
					ty = ((CoordData*)((PathVertex*)terminals[i])->data)->y;
					tz = ((CoordData*)((PathVertex*)terminals[i])->data)->z;
					
					if ( (tx==mx)&&(ty==my)&&(tz==mz))
						found = 1;
				}
				
				/* remove it if we can't find it */
				if (!found) {
					MstVertex	*junk;
					ListElmt *e;
					
					isSteiner = 0;		/* not done yet */
					
					list_rem_next(&T, prev, (void**)(&junk));
					
					/*reset leaves */
					/* initialize */
					for ( e = list_head(&T); e != NULL; e = list_next(e))
					{
						mst_vertex = list_data(e);
						mst_vertex->is_leaf = 1;
					}
					
					/* for each node, set the parent is_leaf to 0.   Then, all leaves will remain */
					for (e = list_head(&T); e != NULL; e = list_next(e))
					{
						mst_vertex = list_data(e);
						
						if (mst_vertex->parent != NULL)
							mst_vertex->parent->is_leaf = 0;
					}
					
					/* start over at beginning of list */
					prev = list_head(&T);
					element = list_next(prev);
					
					
				}
				else {
					prev = list_next(prev);
					element = list_next(element);
				}
			}
			else {
				prev = list_next(prev);
				element = list_next(element);
			}
		}
	}



	/* we can further eliminate vias that connect to a terminal leaf.   These are not neccessary*/
	isSteiner = 0;
	while (!isSteiner) {
		ListElmt	*prev;
		
		/* assume we have it */
		isSteiner = 1;
		
		/* check if each leaf is a terminal */
		prev = list_head(&T);
		element = list_next(prev);
		while (element != NULL) {
			short   int mx,my,mz;
			
			mst_vertex = list_data(element);
			mx = ((CoordData*)((MstVertex*)mst_vertex)->data)->x;
			my = ((CoordData*)((MstVertex*)mst_vertex)->data)->y;		
			mz = ((CoordData*)((MstVertex*)mst_vertex)->data)->z;					
			
			if (mst_vertex->is_leaf) {
				int remove;
				short   int	px,py;
				
				remove = 0;
				px = ((CoordData*)((MstVertex*)mst_vertex->parent)->data)->x;
				py = ((CoordData*)((MstVertex*)mst_vertex->parent)->data)->y;				
				
				if ((px == mx)&&(py == my))
					remove = 1;
				
				
				/* remove it if neccessary */
				if (remove) {
					MstVertex	*junk;
					ListElmt *e;
					
					isSteiner = 0;		/* not done yet */
					
					list_rem_next(&T, prev, (void**)(&junk));
					
					/*reset leaves */
					/* initialize */
					for ( e = list_head(&T); e != NULL; e = list_next(e))
					{
						mst_vertex = list_data(e);
						mst_vertex->is_leaf = 1;
					}
					
					/* for each node, set the parent is_leaf to 0.   Then, all leaves will remain */
					for (e = list_head(&T); e != NULL; e = list_next(e))
					{
						mst_vertex = list_data(e);
						
						if (mst_vertex->parent != NULL)
							mst_vertex->parent->is_leaf = 0;
					}
					
					/* start over at beginning of list */
					prev = list_head(&T);
					element = list_next(prev);
					
					
				}
				else {
					prev = list_next(prev);
					element = list_next(element);
				}
			}
			else {
				prev = list_next(prev);
				element = list_next(element);
			}
		}
	}

	

	*edge_count_OUT = list_size(&T)-1;
	if ((*(SteinerTree_OUT) = (int*) malloc(sizeof(int)*(list_size(&T)-1)))==NULL) {
		printf("gsFH.h : SteinerTree_OUT mem allocation error\n");
		fflush(stdout);
		exit(1);
	}


	i = 0;
	for ( element = list_head(&T); element != NULL; element = list_next(element))
	{
		short   int vx,vy,vz,px,py,pz;
		short   int tx,ty,tz;
		AdjList *a;
		ListElmt *e;
		int		edge_index;
		double	edge_weight;
		
		mst_vertex = list_data(element);
		
		vx = ((CoordData*)mst_vertex->data)->x;
		vy = ((CoordData*)mst_vertex->data)->y;
		vz = ((CoordData*)mst_vertex->data)->z;
		if (mst_vertex->parent != NULL) {
			px = ((CoordData*)mst_vertex->parent->data)->x;
			py = ((CoordData*)mst_vertex->parent->data)->y;		
			pz = ((CoordData*)mst_vertex->parent->data)->z;		
			edge_weight = mst_vertex->weight;
		}
		else {
			px = -1;
			py = -1;
			pz = -1;
		}
		

			
		a = (AdjList*)Edge2Adj[vx][vy][vz];
		
		tx = ((CoordData*)((PathVertex*)(a->vertex))->data)->x;
		ty = ((CoordData*)((PathVertex*)(a->vertex))->data)->y;
		tz = ((CoordData*)((PathVertex*)(a->vertex))->data)->z;
		
		for ( e = list_head(&(a->adjacent)); e != NULL; e = list_next(e) ) {
			PathVertex *p;
		
			p = (PathVertex*)list_data(e);
			
			tx = ((CoordData*)p->data)->x;
			ty = ((CoordData*)p->data)->y;
			tz = ((CoordData*)p->data)->z;
			
			if ((tx == px)&&(ty == py)&&(tz == pz)) {	/*found*/
				edge_index = p->index;
				(*SteinerTree_OUT)[i] = edge_index;
				i++;
			}
		}
	}


	

	/*--------------------------------------------------------------------------------------*/
	/*			Clean Up																	*/
	/*--------------------------------------------------------------------------------------*/
	
	


	/* destroy distance network*/
	graph_destroy(&ND);
	
	/* deystroy all the shortest path lists*/
	for (i = 0; i < NO_TERMINALS; i++) 
	   for (j = 0; j < NO_TERMINALS; j++)
		   if ( i != j ) 
			   list_destroy(&sPaths[i][j]);
	
	/* destroy the pointers to the shortest path lists*/
	for (i = 0; i < NO_TERMINALS; i++) 
		free(sPaths[i]);
	free(sPaths);

	/* destroy the minimum spanning tree*/
	list_destroy(&TD);

	/* destroy grid spanning tree*/
	graph_destroy(&NTD);

	/* destroy the terminal list*/
	for (i = 0; i < NO_TERMINALS; i++) {
		path_vertex_free(terminals[i]);
	}
	free(terminals);
	
	/*destroy the steiner tree*/
	list_destroy(&T);

	return 0;
}
#endif
