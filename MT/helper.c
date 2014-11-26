/*--------------------------------------------------------------------------------------*/
/* helper.c																				*/
/*	Implements routines that are used in many areas of code.   This is pretty bad and	*/
/*	most of these things could definetely be placed elsewhere.							*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/



#ifndef HELPER_C_
#define HELPER_C_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "list.h"
#include "set.h"
#include "graph.h"
#include "graphalg.h"
#include "helper.h"

/*--------------------------------------------------------------------------------------*/
/*						HELPER FUNCTIONS												*/
/*--------------------------------------------------------------------------------------*/

void *path_vertex_free(PathVertex *p) {
	free(p->data);
	free(p);
	return(NULL);
}

void *mst_vertex_free(MstVertex *m) {
	free(m->data);
	free(m);
	return(NULL);
}


void *sPath_vertex_free(sPathData *s) {
	free(s->vertex);
	free(s->parent);
	free(s);
	return(NULL);
}



/*	used to match verticies in the graph routines.   A pointer to this function is passed
	when we call graph_init																*/
int match_coord( const void *grid1, const void *grid2) {	
	int x1,x2,y1,y2,z1,z2;
	
	  x1 = ((CoordData*)(((PathVertex*)grid1)->data))->x;
	  x2 = ((CoordData*)(((PathVertex*)grid2)->data))->x;
	  
	  y1 = ((CoordData*)(((PathVertex*)grid1)->data))->y;
	  y2 = ((CoordData*)(((PathVertex*)grid2)->data))->y;
	  
	  z1 = ((CoordData*)(((PathVertex*)grid1)->data))->z;
	  z2 = ((CoordData*)(((PathVertex*)grid2)->data))->z;

	  
	  return ((x1 == x2)&&(y1 == y2)&&(z1 == z2));
}

void my_itoa(int n, char *n_string) {
	
	int i,digits,tmp;
	
	digits = (int)log10(n) + 1;
	
	if (n == 0) {
		tmp = 0;
		n_string[0] = '0';
		n_string[1] = '\0';
	}
	else {
		int count;
		
		count = 0;
		
		for (i = digits - 1; i >= 0; i--) {
			tmp = n / (int)pow(10,i);
			n_string[count] = tmp + '0';
			count++;
			n -= tmp * (int)pow(10,i);
		}
		n_string[count] = '\0';
	}
}



/* search for an edge in the graph.   Return the cost of the edge */
double find_edge_weight( Graph *graph, CoordData *u, CoordData *v ) {
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
        return 0;
    }
    else {
        
        /* find vertex v in u's adjacency list */
        for (e = list_head(&adjlistPtr->adjacent); e != NULL; e = list_next(e)) {
            PathVertex *tmp;
            int tx,ty,tz;
            
            tmp = list_data(e);
            
            tx = ((CoordData*)((PathVertex*)tmp)->data)->x;
            ty = ((CoordData*)((PathVertex*)tmp)->data)->y;
            tz = ((CoordData*)((PathVertex*)tmp)->data)->z;
            
            if ((tx == v->x)&&(ty == v->y)&&(tz == v->z)) 
                return ((PathVertex*)tmp)->weight;
        }
             
        printf("Couldn't locate v in u's adjacency list\n");
        return 0;
    }
}

int set_edge_weight( Graph *graph, CoordData *u, CoordData *v, int w) {

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
            int tx,ty,tz;
            
            tmp = list_data(e);
            
            tx = ((CoordData*)((PathVertex*)tmp)->data)->x;
            ty = ((CoordData*)((PathVertex*)tmp)->data)->y;
			tz = ((CoordData*)((PathVertex*)tmp)->data)->z;
            
            if ((tx == v->x)&&(ty == v->y)&&(tz == v->z)) {
                tmp->weight = w;
                return 0;
            }
        }
             
        printf("Couldn't locate v (%d,%d) in u's (%d,%d) adjacency list\n",u->x,u->y,v->x,v->y);
        return -1;
    }
}
            
            
    
    

/* each time we call shortest(...) we must copy out the shortest path data.   This is because
shortest() stores the shortest path data within the graph itself.   So, all paths holds is a list of
pointers to verticies in the graph.   It is within the graph structure that the shortest path data is 
contained, so if we call shortest more than once, the graph is re-initialized and all the previous
shortest path computation is lost */
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

void print_weights(Graph *g) {
	ListElmt	*e;
	PathVertex	*v, *a;
	int			vx,vy,vz,ax,ay,az;
	
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		ListElmt	*adj;
		
		v = (PathVertex*)((AdjList*)list_data(e))->vertex;
		
		vx = ((CoordData*)v->data)->x;
		vy = ((CoordData*)v->data)->y;
		vz = ((CoordData*)v->data)->z;
		
		printf("[%d,%d,%d] : ",vx,vy,vz);
		
		for (adj = list_head(&((AdjList*)list_data(e))->adjacent); adj != NULL; adj = list_next(adj)) {
			double	weight;
			
			a = (PathVertex*)list_data(adj);
			weight = a->weight;
			
			ax = ((CoordData*)a->data)->x;
			ay = ((CoordData*)a->data)->y;
			az = ((CoordData*)a->data)->z;
			
			printf("%d:(%d,%d,%d)[%.2f]->",a->index,ax,ay,az,weight);
		}
		
		printf("NULL\n");
	}
	
		
}		


void edge_Histogram(Graph *G, double *f, int base, int size, const char *prefix, int extension) {
	ListElmt	*e;
	PathVertex	*v, *a;
	int			vx,vy,vz,ax,ay,az;
	char		filename[64];
	char		numberstring[10];
	FILE		*fp;
	int			i;
	
	strcpy(filename,prefix);
	strcat(filename,"_");
	my_itoa(extension,numberstring);
	strcat(filename,numberstring);
	strcat(filename,".txt");
	
	if ((fp = fopen( filename, "w")) == NULL){
		printf("can't open file\n");
		exit(1);
	}
	
	/* for each edge */
	for (i = base; i <= size; i++) {
		int found;
		
		/* find edge i in adjacency list */
		e = list_head(&G->adjlists);
		found = 0;
		while ((e != NULL)&&(!found)) {
			ListElmt	*adj;
			
			v = (PathVertex*)((AdjList*)list_data(e))->vertex;
			
			vx = ((CoordData*)v->data)->x;
			vy = ((CoordData*)v->data)->y;
			vz = ((CoordData*)v->data)->z;
			
			adj = list_head(&((AdjList*)list_data(e))->adjacent);
			while( (adj != NULL)&&(!found) ) {
				
				a = (PathVertex*)list_data(adj);
				
				if (a->index == i) {			/* found */

					found = 1;
					
					ax = ((CoordData*)a->data)->x;
					ay = ((CoordData*)a->data)->y;
					az = ((CoordData*)a->data)->z;
					
					fprintf(fp,"%d %d %E\n",vx+ax,vy+ay,f[i]);
				}
				else 
					adj = list_next(adj);
			}
			
			e = list_next(e);
		}
		if ((e == NULL)&&(!found)) 
			printf("Can't find edge %d.\nCheck arguments to function edge_Histogram.\n",i);
	}
	
	fclose(fp);
}




int count_vias(List *L) {
	ListElmt	*element;
	int			vx,vy,vz,px,py,pz;
	MstVertex	*mst_vertex;
	int			vias;
	
	vias = 0;
	for ( element = list_head(L); element != NULL; element = list_next(element)) {
	
		mst_vertex = (MstVertex*) list_data(element);
		
		vx = ((CoordData*)mst_vertex->data)->x;
		vy = ((CoordData*)mst_vertex->data)->y;
		vz = ((CoordData*)mst_vertex->data)->z;
		
		if (mst_vertex->parent != NULL) {
				px = ((CoordData*)mst_vertex->parent->data)->x;
				py = ((CoordData*)mst_vertex->parent->data)->y;		
				pz = ((CoordData*)mst_vertex->parent->data)->z;
				
				if ((vx == px)&&(vy == py))
					vias++;
		}
	}
	
	return(vias);
}	

int tree_length(List *L) {
	ListElmt	*element;
	int			vx,vy,vz,px,py,pz;
	MstVertex	*mst_vertex;
	int			length;
	
	length = 0;
	for ( element = list_head(L); element != NULL; element = list_next(element)) {
		mst_vertex = (MstVertex*) list_data(element);
		vx = ((CoordData*)mst_vertex->data)->x;
		vy = ((CoordData*)mst_vertex->data)->y;
		vz = ((CoordData*)mst_vertex->data)->z;
		
		if (mst_vertex->parent != NULL) {
			px = ((CoordData*)mst_vertex->parent->data)->x;
			py = ((CoordData*)mst_vertex->parent->data)->y;		
			pz = ((CoordData*)mst_vertex->parent->data)->z;
			
			if (pz == vz)
				length++;
		}
	}
	
	return(length);
}

void sync_weights(Graph *g, double *w) {
	ListElmt	*e;
	PathVertex	*a,*p;
	
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		ListElmt	*adj;
		int px,py,pz;
		
		p = (PathVertex*)((AdjList*)list_data(e))->vertex;
		
		px = ((CoordData*)p->data)->x;
		py = ((CoordData*)p->data)->y;
		pz = ((CoordData*)p->data)->z;
		
		for (adj = list_head(&((AdjList*)list_data(e))->adjacent); adj != NULL; adj = list_next(adj)) {
			int	index;
			int ax,ay,az;
			
			
			a = (PathVertex*)list_data(adj);
			index = a->index;
			
			ax = ((CoordData*)a->data)->x;
			ay = ((CoordData*)a->data)->y;
			az = ((CoordData*)a->data)->z;
						
			a->weight = w[index];
		}
	}
	
		
}		

void set_Edge2Adj(Graph *g, AdjList****a) {
	ListElmt		*e;
	PathVertex		*p;
	int				px,py,pz;
	
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		
		p = (PathVertex*)((AdjList*)list_data(e))->vertex;
		
		px = ((CoordData*)p->data)->x;
		py = ((CoordData*)p->data)->y;
		pz = ((CoordData*)p->data)->z;

		a[px][py][pz] = list_data(e);
	}
}
		
void set_edge_lookup(Graph *g, EdgeData *L) {
	ListElmt		*e;
	PathVertex		*p,*a;
	CoordData		*pc,*ac;
	int				px,py,pz,ax,ay,az;
		
	for (e = list_head(&g->adjlists); e != NULL; e = list_next(e)) {
		ListElmt	*adj;
		
		p = (PathVertex*)((AdjList*)list_data(e))->vertex;

		px = ((CoordData*)p->data)->x;
		py = ((CoordData*)p->data)->y;
		pz = ((CoordData*)p->data)->z;
		
		for (adj = list_head(&((AdjList*)list_data(e))->adjacent); adj != NULL; adj = list_next(adj)) {
			int	index;
			
			a = (PathVertex*)list_data(adj);
			index = a->index;
			
			ax = ((CoordData*)a->data)->x;
			ay = ((CoordData*)a->data)->y;
			az = ((CoordData*)a->data)->z;
			
			pc = (CoordData*) malloc(sizeof(CoordData));
			ac = (CoordData*) malloc(sizeof(CoordData));
			
			if ((pc==NULL)||(ac==NULL)) {
				printf("helper.c : set edge lookup mem allocation error\n");
				fflush(stdout);
				exit(1);
			}
			
			pc->x = px;
			pc->y = py;
			pc->z = pz;
			
			ac->x = ax;
			ac->y = ay;
			ac->z = az;
			
			L[index].u = pc;
			L[index].v = ac;
			
		}
	}
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
			
			if ( ((dx==1)||(dx==-1))&&(dy==0) ) /* horizontal edge*/
				(*c)[index] = hc;
			if ( ((dy==1)||(dy==-1))&&(dx==0) ) /* vertical edge*/
				(*c)[index] = vc;
		}
	}
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
			
			if ( ((dx==1)||(dx==-1))&&(dy==0) ) /* horizontal edge*/
				(*c)[index] = hc;
			if ( ((dy==1)||(dy==-1))&&(dx==0) ) /* vertical edge*/
				(*c)[index] = vc;
		}
	}
}

#endif
