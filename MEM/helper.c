/*--------------------------------------------------------------------------------------*/
/* helper.c																				*/
/*	Implements routines that are used in many areas of code.   This is pretty bad and	*/
/*	most of these things could definetely be placed elsewhere.							*/
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


/*--------------------------------------------------------------*/
/*	Purpose: frees path vertices								*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(PathVertex*) p : vertex to be free						*/
/*--------------------------------------------------------------*/
void *path_vertex_free(PathVertex *p) {
	free(p->data);
	free(p);
	return(NULL);
}


/*--------------------------------------------------------------*/
/*	Purpose: frees MST vertices									*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(MSTVertex*) m : vertex to be free						*/
/*--------------------------------------------------------------*/
void *mst_vertex_free(MstVertex *m) {
	free(m->data);
	free(m);
	return(NULL);
}

/*--------------------------------------------------------------*/
/*	Purpose: used to match vertices in other routines			*/
/*  Returns: 1 if grid1 == grid2, 0 otherwise					*/
/*  Args   :												    */
/*		(PathVertex*) grid1,grid2 : vertices to be compared		*/
/*--------------------------------------------------------------*/
int match_coord( const void *grid1, const void *grid2 ) {
	return( (((CoordData*)(((PathVertex*)grid1)->data))->x == ((CoordData*)(((PathVertex*)grid2)->data))->x) &&
			(((CoordData*)(((PathVertex*)grid1)->data))->y == ((CoordData*)(((PathVertex*)grid2)->data))->y) &&
			(((CoordData*)(((PathVertex*)grid1)->data))->z == ((CoordData*)(((PathVertex*)grid2)->data))->z) );
}

int match_coord2( const void *grid1, const void *grid2) {
	return( (((CoordData2*)(((PathVertex*)grid1)->data))->x == ((CoordData2*)(((PathVertex*)grid2)->data))->x) &&
			(((CoordData2*)(((PathVertex*)grid1)->data))->y == ((CoordData2*)(((PathVertex*)grid2)->data))->y) );
}


/*--------------------------------------------------------------*/
/*	Purpose: converts an integer to a string					*/
/*  Returns: (char*) n_string									*/
/*  Args   :												    */
/*		(int) n : integer to be stringified						*/
/*		(char*) n_string : string version of n					*/
/*--------------------------------------------------------------*/
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


#endif
