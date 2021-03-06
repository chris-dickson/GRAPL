static void *path_vertex_free(PathVertex *p) {
	free(p->data);
	free(p);
	return(NULL);
}

static int qsort_match(const void *a1, const void *a2) {
	if (*(int*)a1 < *(int*)a2)
		return -1;
	else if (*(int*)a1 > *(int*)a2)
		return 1;
	return 0;
}


static int match_coord( const void *grid1, const void *grid2) {	
	int x1,x2,y1,y2,z1,z2;
	
	x1 = ((CoordData*)(((PathVertex*)grid1)->data))->x;
	x2 = ((CoordData*)(((PathVertex*)grid2)->data))->x;
	
	y1 = ((CoordData*)(((PathVertex*)grid1)->data))->y;
	y2 = ((CoordData*)(((PathVertex*)grid2)->data))->y;
	
	z1 = ((CoordData*)(((PathVertex*)grid1)->data))->z;
	z2 = ((CoordData*)(((PathVertex*)grid2)->data))->z;
	
	
	return ((x1 == x2)&&(y1 == y2)&&(z1 == z2));
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// gen_gridgraph( ... )
//		This function generates a grid graph with vertices [0..GRID_WIDTH) x [0..GRID_HEIGHT).
//	Additionally, it adds 4 edges to each internal vertex. 2 edges to corner vertices, and 3 edges to
//	each vertex on the boundary.   All edges added are given weight 1.   
//	The graph is returned as the first parameter to the function.   There is no need to allocate/init
//	a graph before calling this function
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int gen_gridgraph(Graph *grid_graph, int GRID_WIDTH, int GRID_HEIGHT) {
	
	PathVertex		*path_vertex, *path_start, *path_vertex1, *path_vertex2;
	CoordData		*coord;
	int				i,j;
	int				index;
	
	printf("Creating Grid Graph\n");
	
	/* initialize a grid graph */
	graph_init(grid_graph, match_coord, (void*)path_vertex_free);
	printf("\tInitialized grid graph..\n");
	
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
            if (graph_ins_vertex(grid_graph, path_vertex) != 0) {
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
                
                if (graph_ins_edge(grid_graph, path_vertex1, path_vertex2) != 0)
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
                
                if (graph_ins_edge(grid_graph, path_vertex1, path_vertex2) != 0)
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
                
                if (graph_ins_edge(grid_graph, path_vertex1, path_vertex2) != 0)
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
                
                if (graph_ins_edge(grid_graph, path_vertex1, path_vertex2) != 0)
                    return 1;
				
				free(path_vertex1->data);
				free(path_vertex1);
				
            }
			
			
        }
    }/* finished inserting edges */
	printf("Finished\n");
	return 0;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// grid2VL( ... )
//		This function converts a grid graph G to a virtual (2 layer) graph H.  There is no need to 
//	initialize/allocate H, but G should already be created.   If we are removing vertices (ie/ if the
//	graph contains holes, this should be done BEFORE calling this function, as in this function, each 
//	edge is enumerated with an integer.   Since we convert G to H, all edges in H are initially weight
//	1.   Use sync_weights(...) to change the weights of edges after calling grid2VL().
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void grid2VL(Graph *G, Graph *H, int width, int height) { 
	ListElmt	*elmt;
	int			index,via;
	
	
	printf("Creating Virtual Layer Graph\n");
	
	graph_init(H, match_coord, (void*)path_vertex_free);
	
	// original edges enumrated from 1, vias enumerated after original edges.
	// also, since our graph is undirected, we have to half edges for each edge in undirected graph, 
	// thus we start vias from ecount/2
	index = 1;
	via = (G->ecount / 2) + 1;
	
	
	// for each vertex in G, there will be two vertices in H.  (One per layer)
	printf("\tInserting vertices and vias into multi-layer graph....");
	fflush(stdout);
	for (elmt = list_head(&G->adjlists); elmt != NULL; elmt = list_next(elmt)) {
		PathVertex	*v0,*v1,*v;		// v = vertex from g, v1,v2 = new vertices to insert in H
		CoordData	*v0c,*v1c;
		
		v = ((AdjList*)list_data(elmt))->vertex;
		
		// insert both verticies
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
		
		// horizontal vertex
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;
		
		// vertical vertex
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		graph_ins_vertex(H,v0);
		graph_ins_vertex(H,v1);
		
		// now we insert the via between v0 and v1
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
		
		// set the weights to 1, temporarily
		v1->weight	 = 1;
		
		// set via index
		v1->index = via;
		
		// horizontal vertex
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;
		
		// vertical vertex
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		
		// forward edge
		graph_ins_edge(H,v0,v1);
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
		
		// set the weights to 1, temporarily
		v0->weight	 = 1;
		
		// via index
		v0->index = via;
		
		// horizontal vertex
		v0c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v0c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v0c->z = 0;
		
		// vertical vertex
		v1c->x = ((CoordData*)((PathVertex*)v)->data)->x;
		v1c->y = ((CoordData*)((PathVertex*)v)->data)->y;
		v1c->z = 1;
		
		// forward edge
		graph_ins_edge(H,v1,v0);
		free(v1c);
		free(v1);
		
		// next via
		via++;
		
	}
	printf("Finished\n");
	fflush(stdout);
	
	
	printf("\tInserting edges into multi-layer graph..");
	fflush(stdout);
	// now we insert edges into the multi layer graph
	for (elmt = list_head(&G->adjlists); elmt != NULL; elmt = list_next(elmt)) {
		PathVertex	*v0,*v1,*v,*next;		// v = vertex from g, v1,v2 = new vertices to insert in H
		CoordData	*v0c,*v1c;
		ListElmt	*e;
		int			vx,vy,vz;
		
		//printf(".");
		//fflush(stdout);
		
		v = ((AdjList*)list_data(elmt))->vertex;
		
		vx = ((CoordData*)((PathVertex*)v)->data)->x;
		vy = ((CoordData*)((PathVertex*)v)->data)->y;
		vz = ((CoordData*)((PathVertex*)v)->data)->z;
		
		// go through each vertex adjacent to v
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
			
			// find out if it's a horizontal or vertical edge
			if (vx == nx) layer = 1;
			else if (vy == ny) layer = 0;
			else assert(0);
			
			// insert forward edge from v0 to v1
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
			
			graph_ins_edge(H,v0,v1);
			
			// free the lookup variable
			free(v0c);
			free(v0);
			
			
			// insert the backedge v1 to v0
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
			
			exists = graph_ins_edge(H,v1,v0);
			
			if (!exists) {
				((PathVertex*)next)->index = index;		// set index in orginal graph G
				index++;
			}
			
			// free the lookup variable
			free(v1c);
			free(v1);
			
		}
	}
	printf("finished\n");
	fflush(stdout);
	
	
}

void set_capacity_lookup(int **c, Graph *g, int hc, int vc) {
	ListElmt		*e;
	PathVertex		*p,*a;
	int				px,py,pz,ax,ay,az;
	
	
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		ListElmt	*adj;
		
		p = (PathVertex*)((AdjList*)list_data(e))->vertex;
		
		px = ((CoordData*)p->data)->x;
		py = ((CoordData*)p->data)->y;
		pz = ((CoordData*)p->data)->z;
		
		for (adj = list_head(&((AdjList*)list_data(e))->adjacent); adj != NULL; adj = list_next(adj)) {
			int	index;
			int dx,dy,dz;
			
			a = (PathVertex*)list_data(adj);
			index = a->index;
			
			ax = ((CoordData*)a->data)->x;
			ay = ((CoordData*)a->data)->y;
			az = ((CoordData*)a->data)->z;
			
			dx = ax-px;
			dy = ay-py;
			dz = az-pz;
			
			if ( ((dx==1)||(dx==-1))&&(dy==0) ) // horizontal edge
				(*c)[index] = hc;
			if ( ((dy==1)||(dy==-1))&&(dx==0) ) // vertical edge
				(*c)[index] = vc;
		}
	}
}

int U(Net *T, int edge_i, int tree_j) {
	int count;
	
	for (count = 0; count < T[tree_j].edge_count; count++)
		if (T[tree_j].SteinerTree[count] == edge_i)
			return 1;
	
	return 0;
}

void set_length_lookup(double **c, Graph *g, double hc, double vc) {
	ListElmt		*e;
	PathVertex		*p,*a;
	int				px,py,pz,ax,ay,az;
	
	
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		ListElmt	*adj;
		
		p = (PathVertex*)((AdjList*)list_data(e))->vertex;
		
		px = ((CoordData*)p->data)->x;
		py = ((CoordData*)p->data)->y;
		pz = ((CoordData*)p->data)->z;
		
		for (adj = list_head(&((AdjList*)list_data(e))->adjacent); adj != NULL; adj = list_next(adj)) {
			int	index;
			int dx,dy,dz;
			
			a = (PathVertex*)list_data(adj);
			index = a->index;
			
			ax = ((CoordData*)a->data)->x;
			ay = ((CoordData*)a->data)->y;
			az = ((CoordData*)a->data)->z;
			
			dx = ax-px;
			dy = ay-py;
			dz = az-pz;
			
			if ( ((dx==1)||(dx==-1))&&(dy==0) ) // horizontal edge
				(*c)[index] = hc;
			if ( ((dy==1)||(dy==-1))&&(dx==0) ) // vertical edge
				(*c)[index] = vc;
		}
	}
}


// use this one only when the list of edge indices is sorted
int U_sorted(Net *T, int edge_i, int tree_j) {
	
	if (bsearch( &edge_i , T[tree_j].SteinerTree , T[tree_j].edge_count , sizeof(int), &qsort_match ) == NULL)
		return 0;
	else
		return 1;
}
	


void print_all_nets(FILE *fp,Net *T, int N, int K) {
	int l;
	
	for (l = 0; l <= N; l++) {
		int k;
		
		fprintf(fp,"Iteration %d\n",l);
		fprintf(fp,"------------\n\n");
		for (k = 1; k <= K; k++) {
			int i;
			
			fprintf(fp,"\tNet %d\n",k);
			fprintf(fp,"\t\tno_terminals=%d\n",T[l*K + k].no_terminals);
			fprintf(fp,"\t\tterminal set:\n");
			for (i = 0; i < T[l*K + k].no_terminals; i++)
				fprintf(fp,"\t\t\t%d : (%d,%d,%d)\n",i,T[l*K + k].terminals[i].x,T[l*K + k].terminals[i].y,T[l*K + k].terminals[i].z);
			fprintf(fp,"\n");
			fprintf(fp,"\t\tedge_count=%d\n",T[l*K + k].edge_count);
			fprintf(fp,"\t\tvias=%d\n",T[l*K + k].vias);
			fprintf(fp,"\t\ttree_cost=%f\n",T[l*K + k].tree_cost);
			fprintf(fp,"\t\tsteiner tree:\n");
			for (i = 0; i < T[l*K + k].edge_count; i++) 
				fprintf(fp,"\t\t\t%d\n",T[l*K + k].SteinerTree[i]);
			fprintf(fp,"\n\n");
			
		}
		fprintf(fp,"\n");
	}
	
}


