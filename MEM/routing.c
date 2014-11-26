/*--------------------------------------------------------------------------------------*/
/* routing.c																			*/
/*	This file implements some standard routines for our routing problem.  It contains	*/
/* debug routines, as well as routines to parse the benchmarks.							*/
/*--------------------------------------------------------------------------------------*/

#ifndef ROUTING_C_
#define ROUTING_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "graphtype.h"
#include "routing.h"
#include "set.h"


/*--------------------------------------------------------------*/
/*	Purpose: prints entire net array to file fp					*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(FILE*) fp  : file to print to							*/
/*		(Net*) T :  our net array								*/
/*		(int) N : max number of iterations						*/
/*		(int) K : number of nets								*/
/*--------------------------------------------------------------*/
void print_all_nets(FILE *fp,Net *T, int N, int K) {
	int l;
	
	for (l = 0; l <= N; l++) {
		int k;
		
		fprintf(fp,"Iteration %d\n",l);
		fprintf(fp,"------------\n\n");
		for (k = 1; k <= K; k++) {
			int i;
			
			fprintf(fp,"\tNet %d\n",k);
			fprintf(fp,"\tID %s\n",T[l*K + k].ID);
			fprintf(fp,"\t\tno_terminals=%d\n",T[l*K + k].no_terminals);
			fprintf(fp,"\t\tterminal set:\n");
			for (i = 0; i < T[l*K + k].no_terminals; i++)
				fprintf(fp,"\t\t\t%d : (%d,%d,%d)\n",i,T[l*K + k].terminals[i].x,T[l*K + k].terminals[i].y,T[l*K + k].terminals[i].z);
			fprintf(fp,"\n");
			
			if (T[l*K + k].routed) {
				fprintf(fp,"\t\tedge_count=%d\n",T[l*K + k].edge_count);
				fprintf(fp,"\t\tvias=%d\n",T[l*K + k].vias);
				fprintf(fp,"\t\ttree_cost=%f\n",T[l*K + k].tree_cost);
				fprintf(fp,"\t\tsteiner tree:\n");
				for (i = 0; i < T[l*K + k].edge_count; i++) 
					fprintf(fp,"\t\t\t%d\n",T[l*K + k].SteinerTree[i]);
				fprintf(fp,"\n\n");
			}
				
		}
		fprintf(fp,"\n");
	}
	
}


/*--------------------------------------------------------------*/
/*	Purpose: prints net[k]  to file fp							*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(FILE*) fp  : file to print to							*/
/*		(Net*) T :  our net array								*/
/*		(int) k : net to print									*/
/*--------------------------------------------------------------*/

void print_one_net(FILE *fp, Net *T, int k) {
	int i;
	
	fprintf(fp,"\tNet %d\n",k);
	fprintf(fp,"\tID %s\n",T[k].ID);
	fprintf(fp,"\t\tno_terminals=%d\n",T[k].no_terminals);
	fprintf(fp,"\t\tterminal set:\n");
	for (i = 0; i < T[k].no_terminals; i++)
		fprintf(fp,"\t\t\t%d : (%d,%d,%d)\n",i,T[k].terminals[i].x,T[k].terminals[i].y,T[k].terminals[i].z);
	fprintf(fp,"\n");
	
//	if (T[k].routed) {
		fprintf(fp,"\t\tedge_count=%d\n",T[k].edge_count);
		fprintf(fp,"\t\tvias=%d\n",T[k].vias);
		fprintf(fp,"\t\ttree_cost=%f\n",T[k].tree_cost);
		fprintf(fp,"\t\tsteiner tree:\n");
		for (i = 0; i < T[k].edge_count; i++) 
			fprintf(fp,"\t\t\t%d\n",T[k].SteinerTree[i]);
		fprintf(fp,"\n\n");
//	}
}


void print_range_net(FILE *fp, Net *T, int lb, int ub){
}


/*--------------------------------------------------------------*/
/*	Purpose: prints entire net array to file fp					*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(FILE*) fp  : file to print to							*/
/*		(Net*) T :  our net array								*/
/*		(int) N : max number of iterations						*/
/*		(int) K : number of nets								*/
/*--------------------------------------------------------------*/
void dump_Array(double *a, const char *prefix, int extension, int size) {
	char filename[64];
	char numberstring[10];
	FILE *fp;
	int i;
	
	
	
	strcpy(filename,prefix);
	strcat(filename,"_");
	my_itoa(extension,numberstring);
	strcat(filename,numberstring);
	strcat(filename,".txt");
	
	if ((fp = fopen( filename, "w")) == NULL){
		printf("can't open file\n");
		exit(1);
	}
	else {
		for (i = 1; i <= size; i++)
			fprintf(fp,"%d %E\n",i,a[i]);
		fclose(fp);
	}
}


/*--------------------------------------------------------------*/
/*	Purpose: determines if edge e(i) is used in tree T(j)		*/
/*  Returns: 1 if edge(i) is in tree(j), 0 otherwise			*/
/*  Args   :												    */
/*		(Net*) T  : our main data array							*/
/*		(int) edge_i : integer index of edge e(i)				*/
/*		(int) tree_j : integer index of tree T(j)				*/
/*--------------------------------------------------------------*/
int U(Net *T, int edge_i, int tree_j) {
	int count;
	
	for (count = 0; count < T[tree_j].edge_count; count++)
		if (T[tree_j].SteinerTree[count] == edge_i)
			return 1;
	
	return 0;
}

typedef struct point2d_ {
	int x,y;
} point2d;


/*--------------------------------------------------------------*/
/*	Purpose: determines if point2d vars are equal				*/
/*  Returns: 1 if a1 == a2, 0 otherwise							*/
/*  Args   :												    */
/*		(point2d*) a1,a2  : points to be compared				*/
/*--------------------------------------------------------------*/
static int match2d(const void *a1, const void *a2) {
	return ((((point2d*)a1)->x == ((point2d*)a2)->x)&&(((point2d*)a1)->y == ((point2d*)a2)->y));
}


/*--------------------------------------------------------------*/
/*	Purpose: parse through a benchmark file to create data		*/
/*			array T												*/
/*  Returns:													*/
/*		(Net*) T : main data array								*/
/*		(int*) gw,gh : width and height of instance				*/
/*		(int*) hCap,vCap : horizontal and vertical capacities	*/
/*		(int*) num_routed : always the same as number in		*/
/*				benchmark file									*/
/*  Args   :												    */
/*		(char*) fName  : filename of benchmark file				*/
/*					.................							*/
/*	Other Notes:												*/
/*		The file must be pruned!!  We used to have pruning		*/
/*	code integrated into this procedure.  But this is now		*/
/*	depricated.													*/
/*		THIS IS FOR IBM_STYLE BENCHMARKS ONLY!!					*/
/*--------------------------------------------------------------*/
int parseBenchmark(char *fName, int *no_nets, Net **T, int N, int *gw, int *gh, int *hCap, int *vCap, int *num_routed) {
	FILE	*fp;
	int		i,l;
	int		K;
	Net		*temp;
	
	double	AvgNets;
	int		tolerance;
	int		routed;
	double	thrown_out;
	int		net_array_size;
	int	 no_term_table[100];
	 
	
	
	
	AvgNets = 0;
	tolerance = 6;
	routed = 0;
	
	/* open the benchmark file*/
	if ( (fp = fopen(fName,"r")) == NULL) {
		printf("Can't find the file "); printf(fName); printf("\n");
		*no_nets = 0;
		*T = NULL;
		*gw = 0;
		*gh = 0;
		return 1;
	}
	
	/* grid dimensions*/
	fscanf(fp,"%*s %d %d",gw,gh);
	
	/* vertical capacity*/
	fscanf(fp,"%*s %*s %d", vCap);
	
	/* horizontal capacity*/
	fscanf(fp,"%*s %*s %d", hCap);
	
	/* num nets*/
	fscanf(fp,"%*s %*s %d",no_nets);
	
	K = *no_nets;
	
	/* (N+1) iterations total, each of K nets (+1 so we can index iteration 0 from T[1]..T[K]*/
	net_array_size = K*(N+1) + 1;
		
	if ((temp = (Net*)malloc(sizeof(Net)*net_array_size))==NULL){
		printf("benchmark, temp : memory allocation error\n"); fflush(stdout);
		exit(1);
	}
	
	/*begin parsing nets*/
	printf("Parsing nets...");
	for (i = 1; i <= *no_nets; i++) {
		Set			pointSet;
		int			numPoints,j;
		ListElmt		*e;
		int			net_number;
		char		ID[ID_LENGTH];
		
		/* initalize a point set*/
		set_init( &pointSet, match2d, free );
		

		fscanf(fp,"%s %d %d",ID,&net_number,&numPoints);
		
		no_term_table[ numPoints ]++;
		
		
		/* throw all points into a set, duplicates will get removed*/
		for (j = 0; j < numPoints; j++) {
			point2d *tmp;
			int		x,y;
	
			fscanf(fp,"%d %d",&x,&y);
			
			if ((tmp = (point2d*) malloc(sizeof(point2d)))==NULL){
				printf("benchmarks.c : tmp mem allocation error\n");
				fflush(stdout);
				exit(1);
			}
			tmp->x = x;	tmp->y = y;
			
			if (set_insert(&pointSet,tmp) == 1) {
				/* already in set, so destroy tmp */
				free(tmp);
			}
		}
		
		
		/* convert point set to an array of CoordData*/
		if ((temp[i].terminals = (CoordData*) malloc(sizeof(CoordData) * set_size(&pointSet))) == NULL) {
			printf("benchmark, temp[i].terminals : memory allocation error\n"); fflush(stdout);
			fflush(stdout);
			exit(1);
		}
		
		strcpy(temp[i].ID, ID);
		temp[i].net_num = net_number;
		temp[i].no_terminals = set_size(&pointSet);
		temp[i].routed = 0;
		temp[i].fixed = 0;
		AvgNets += set_size(&pointSet);
		if ((set_size(&pointSet) <= tolerance))
			routed++;
		j = 0;
		for (e = list_head(&pointSet); e != NULL; e = list_next(e)) {
			point2d		*data;
			
			data = (point2d*) list_data(e);
			temp[i].terminals[j].x = data->x; temp[i].terminals[j].y = data->y;
			temp[i].terminals[j].z = 0;
			j++;
		}		
		set_destroy(&pointSet);
	}
	printf("done\n");
	
	AvgNets /= *no_nets;
	thrown_out = routed/(*no_nets);
	
	
	*num_routed = routed;
	
	
	for (l = 1; l <= N; l++) 
		memcpy( temp + (l*K) + 1 , temp + 1 , sizeof(Net)*K);
	
	*T = temp;
	
	fclose(fp);

	return 0;
}	



/*--------------------------------------------------------------*/
/*	Purpose: parse through a benchmark file to create data		*/
/*			array T												*/
/*  Returns:													*/
/*		(Net*) T : main data array								*/
/*		(int*) gw,gh : width and height of instance				*/
/*		(int*) hCap,vCap : horizontal and vertical capacities	*/
/*		(double*) hLength, vLength : horizontal and vertical	*/
/*				edge lengths.  (physical lengths)				*/
/*  Args   :												    */
/*		(char*) fName  : filename of benchmark file				*/
/*					.................							*/
/*	Other Notes:												*/
/*		The file must be pruned!!  We used to have pruning		*/
/*	code integrated into this procedure.  But this is now		*/
/*	depricated.													*/
/*		THIS IS FOR YAL_STYLE BENCHMARKS ONLY!!					*/
/*--------------------------------------------------------------*/
void parse_yal_benchmark(char *fName, int *no_nets, Net **T, int N, int *gw, int *gh, int *hCap, int *vCap, double *hLength, double *vLength) {
	FILE *fp;
	Net	 *temp;
	int	 net_array_size,k,l,K;
	int	 no_term_table[100];
	
	fp = fopen(fName,"r");
	
	
	fscanf(fp,"%*s %*s %*s %d",no_nets);	
	K = *no_nets;
			
	fscanf(fp,"%*s %*s %*s %*d");					
	fscanf(fp,"%*s %*s %*s %d",gw);					
	fscanf(fp,"%*s %*s %*s %d",gh);					
	fscanf(fp,"%*s %*s %*s %d",hCap);			
	*vCap = *hCap;
	
	fscanf(fp,"%*s %*s %*s %lf",hLength);
	fscanf(fp,"%*s %*s %*s %lf",vLength);
		
	net_array_size = K*(N+1) + 1;	
	if ((temp = (Net*)malloc(sizeof(Net)*net_array_size))==NULL){
		printf("benchmark, temp : memory allocation error\n"); fflush(stdout);
		exit(1);
	}
	
	for (k = 0; k < 100; k++) 
		no_term_table[k] = 0;
	
	
	for (k = 1; k <= K; k++) {
		int i;
		
		fscanf(fp,"%d, %d:,",&(temp[k].net_num),&(temp[k].no_terminals));
		
		no_term_table[ temp[k].no_terminals ]++;
		
		strcpy(temp[k].ID,"netXXXXX");
		temp[k].terminals = (CoordData*) malloc(sizeof(CoordData)*temp[k].no_terminals);
		temp[k].routed = 0;
		temp[k].fixed = 0;
		
		for (i = 0; i < temp[k].no_terminals; i++) {
			char a_string[16];
			
			fscanf(fp,"%s",a_string);
			sscanf(a_string,"%hd,%hd,",&(temp[k].terminals[i].x),&(temp[k].terminals[i].y));
			temp[k].terminals[i].z = 0;
			
			if (i == temp[k].no_terminals - 1) 
				fgets(a_string,16,fp);			/* pick up the trash*/
			
		}

	}
	fclose(fp);
	
	for (l = 1; l <= N; l++) 
		memcpy( temp + (l*(*no_nets)) + 1 , temp + 1 , sizeof(Net)*(*no_nets));
	

	
	*T = temp;
}


/*--------------------------------------------------------------*/
/*	Purpose: determine which coordinates are a terminal in some	*/
/*			net.												*/
/*  Returns: 													*/
/*		(CoordData**) unique_terminals : pointer to array of	*/
/*			terminals that are used in some net					*/
/*		(int*) no_unique_terminals : size of unique_terminals	*/
/*  Args   :												    */
/*		(Net*) T  : main data array returned from parse_...		*/
/*					.................							*/
/*--------------------------------------------------------------*/
void find_unique_terminals(Net *T, int no_nets, CoordData **unique_terminals, int *no_unique_terminals){
	int			k,i;
	Set			pointSet;
	ListElmt	*e;
	CoordData	*tmp_UT;
	
	set_init( &pointSet, match2d, free );
	
	for (k = 1; k <= no_nets; k++) {
		
		for (i = 0; i < T[k].no_terminals; i++) {
			point2d	*p;
		
			if ((p = (point2d*) malloc(sizeof(point2d)))==NULL) {
				printf("benchmarks.c : tmp mem allocation error\n");
				fflush(stdout);
				exit(1);
			}
			
			p->x = (int)T[k].terminals[i].x;
			p->y = (int)T[k].terminals[i].y;
			
			if (set_insert(&pointSet,p) == 1) {
				/* already in set, so destroy p */
				free(p);
			}
		}
	}
	
		
	*no_unique_terminals = set_size(&pointSet);
	
	if ((tmp_UT = (CoordData*) malloc(sizeof(CoordData)* set_size(&pointSet)))==NULL) {
		printf("benchmark, tmp_UT : memory allocation error\n"); fflush(stdout);
		fflush(stdout);
		exit(1);
	}

	
	i = 0;
	for (e = list_head(&pointSet); e != NULL; e = list_next(e), i++) {
		point2d	*data;
		
		data = (point2d*) list_data(e);
				
		tmp_UT[i].x = (short unsigned int)data->x;
		tmp_UT[i].y = (short unsigned int)data->y;
		tmp_UT[i].z = 0;
	}
	
	list_destroy(&pointSet);
	
	*unique_terminals = tmp_UT;
}


/*--------------------------------------------------------------*/
/*	Purpose: debug routine.   Used to see how pruning affects	*/
/*		number of nets included in an instance					*/
/*  Returns: nothing											*/
/*  Args   :												    */
/*		(Net*) T  : main data array								*/
/*		(int) K : number of nets								*/
/*	Other Notes:												*/
/*		This can only be entered by turning on a #define flag	*/
/*	in the main c file.											*/
/*--------------------------------------------------------------*/
void query_mode(Net *T,int K) {
	char fn[64];
	FILE *fp;
	
	int i,done,tol,count;
	double ratio;
	
	
	printf("Query Mode\n");
	printf("----------\n");
	
	done = 0;
	fp = NULL;
	while (!done) {
		printf("Enter the benchmark file to be read: "); 
		scanf("%s",fn);
		
		if (strcmp(fn,"ls") == 0) 
			system("ls");
		else if (strcmp(fn,"quit") == 0)
			exit(0);
		else 
			if (( fp = fopen(fn,"r") ) == NULL) 
				printf("Can't find the file %s\n",fn);
		else 
			done = 1;
		
	}
	
	while (1) {
		printf("Enter a tolerance (integeral, 0 to quit) : ");
		scanf("%d",&tol);
		
		if (tol == 0) {
			fclose(fp);
			return;
		}
		
		count = 0;
		for (i = 1; i <= K; i++) 
			if (T[i].no_terminals <= tol)
				count++;
		
		printf("\t# nets include    : %d\n",count);
		printf("\t# nets thrown out : %d\n",K-count);
		printf("\t       TOTAL NETS : %d\n",K);
		ratio = (double) (count / K)*100;
		printf("\tTolerance of %d throws out %d percent of the total nets\n\n",tol, (int)(((double)(K-count)/(double)K)*100) );
	}
	
	
}


/*--------------------------------------------------------------*/
/*	Purpose: dump a file containing edge congestion to be read	*/
/*		by visualization code OPTRouteViz	(MATLAB)			*/
/*  Returns: nothing, writes a file determined by prefix		*/
/*		appended by the iteration number						*/
/*  Args   :												    */
/*		(Graph*) G  : pointer to the graph						*/
/*		(double*) f : scaled edge congetsion					*/
/*		(int) base : index of starting constraint (usually 0)	*/
/*		(int) size : number of constraints to include			*/
/*		(const char*) prefix : prefix of file name				*/
/*			ie/ Congestion___.txt								*/
/*		(int) extension : iteration number to be appended to	*/
/*			the filename										*/
/*	Other Notes:												*/
/*		This program computes the dual of the graph.  That is,	*/
/*	edges become vertices in the histogram.						*/
/*--------------------------------------------------------------*/
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
			
			vx = (int)((CoordData*)v->data)->x;
			vy = (int)((CoordData*)v->data)->y;
			vz = (int)((CoordData*)v->data)->z;
			
			adj = list_head(&((AdjList*)list_data(e))->adjacent);
			while( (adj != NULL)&&(!found) ) {
				
				a = (PathVertex*)list_data(adj);
				
				if (a->index == i) {			/* found */
					
					found = 1;
					
					ax = (int)((CoordData*)a->data)->x;
					ay = (int)((CoordData*)a->data)->y;
					az = (int)((CoordData*)a->data)->z;
					
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

void fix_net(Net *T, int i) {
	T[i].fixed = 1;
} 

void
get_next_tree(FILE *fp, Graph *g, double hl, double vl, int **SteinerTree, int *edge_count_OUT, double *tree_cost_OUT, char *ID_OUT) 
{
	char	ID[ID_LENGTH];
	char	a_line[LINELENGTH];
	int		no_edges;
	double	tree_cost;
	int		Edges[MAX_NO_EDGES];
	int		i;
	
	
	no_edges = 0;
	tree_cost = 0.0;
	
	fscanf(fp,"%*s %*s %s",ID);		/* net name: ID				*/
    strcpy(ID_OUT,ID);
	fscanf(fp,"%*s \n");			/* edges:					*/

	fgets(a_line,LINELENGTH,fp);
	while ((strcmp(a_line,"###################################\n"))&&(strcmp(a_line,"\n")))
	{
		char		direction[LINELENGTH];
		int			track,edge;
		CoordData	u,v;
		
		sscanf(a_line,"%s %*s %*s %*s %d %*s %*s %d\n",direction,&track,&edge);
		
		if (!strcmp(direction,"horizontal")) {
			u.x = edge;
			u.y = track;
			u.z = 0;
			
			v.x = edge+1;
			v.y = track;
			v.z = 0;
			
			tree_cost += hl;
			
		}
		else {
			u.x = track;
			u.y = edge;
			u.z = 1;

			v.x = track;
			v.y = edge+1;
			v.z = 1;
			
			tree_cost += vl;
			
		}
		
		/* lookup edge */
		Edges[no_edges] = find_edge_index(g,&u,&v);
		no_edges++;		
		fgets(a_line,LINELENGTH,fp);
	}
		   
    *(SteinerTree) = (int*) malloc(sizeof(int)*no_edges);
    for (i = 0; i < no_edges; i++)
		(*SteinerTree)[i] = Edges[i];   


	(*edge_count_OUT) = no_edges;
    (*tree_cost_OUT) = tree_cost;
		   
}
			

			
			   
		
		
		
	
	


#endif
