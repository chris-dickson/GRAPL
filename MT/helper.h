/*--------------------------------------------------------------------------------------*/
/* helper.h																				*/
/*	Header for auxilliary routines that are used in many areas of code.					*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/



#ifndef HELPER_H_
#define HELPER_H_

#define MAX_TERMINALS 20
#define MAX_WIDTH 256
#define MAX_HEIGHT 64

#include "list.h"
#include "graph.h"
#include "graphalg.h"

/* coordinate data in graphs */
typedef struct CoordData_ {
		
	int					x,
						y,
						z;
	
} CoordData;

/* edge data structure */
typedef struct EdgeData_ {

CoordData              *u,*v;
double					weight;

} EdgeData;

typedef struct sPathData_ {
	CoordData			*vertex,
						*parent;
	
	double				d;
} sPathData;


void *path_vertex_free(PathVertex *p);
void *mst_vertex_free(MstVertex *m);
void *sPath_vertex_free(sPathData *s);
int match_coord( const void *grid1, const void *grid2);
void my_itoa(int n, char *n_string);

double find_edge_weight( Graph *graph, CoordData *u, CoordData *v );
int set_edge_weight( Graph *graph, CoordData *u, CoordData *v, int w);
int copy_sPath(List *path, List *sPath, CoordData *start);

void print_weights(Graph *g);
void edge_Histogram(Graph *G, double *f, int base, int size, const char *prefix, int extension);
int count_vias(List *L);
int tree_length(List *L);
void sync_weights(Graph *g, double *w);
void set_Edge2Adj(Graph *g, AdjList****a);
void set_edge_lookup(Graph *g, EdgeData *L);
void set_capacity_lookup(int **c, Graph *g, int hc, int vc);
void set_length_lookup(double **c, Graph *g, double hc, double vc);

#endif
