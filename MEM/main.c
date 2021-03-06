/*--------------------------------------------------------------------------------------*/
/* main.c																				*/
/*	This is our main function.   The purpose of this code is to approximate the			*/
/* solution to an ILP formulation of the global routing problem.   For more infor		*/
/* on this algorithm, please see the corresponding paper.								*/
/*																						*/
/*	Use:																				*/
/*		route inputfile.txt args.txt													*/
/*																						*/
/*	Where inputfile.txt is a text file that is either a Labyrinth style benchmark or an */
/* MCNC style benchmark.																*/
/*																						*/
/* The argument is a list of parameters settings, one setting per line, as follows:		*/
/*		setting=value;																	*/
/*																						*/
/*	Where setting is the parameter we wish to set, and value is the value which it is	*/
/* to be set to.   We do not allow spaces between the equal sign or the semi colon.		*/
/*																						*/
/*	Parameters:																			*/
/*																						*/
/*		N				:	Maximum number of iterations (default = 10)					*/
/*		EPS				:	The value of epsilon (default = 0.001)						*/
/*		TAU_HEURISTIC	:	The way that our step length is calculated					*/
/*							1 : default, as in paper									*/
/*							2 :															*/
/*							3 :															*/
/*							4 :															*/
/*							5 :															*/
/*							6 :															*/
/*							7 :															*/
/*							8 :															*/
/*							9 :															*/
/*							10:															*/
/*							11:															*/
/*		ALPHA			:	The artificial weight of the original grid edges.   A value */
/*							for this parameter puts a higher importance to minimizing	*/
/*							total wirelength. (default = 0.5)							*/
/*		BETA			:	The artificial weight of the vias in the VL graph.  A high	*/
/*							value for beta should reduce the number of vias in the		*/
/*							solution.	(default = 0.5)									*/
/*		G				:	Our guessed objective value. 								*/
/*							Default value for this is the maximum number possible		*/
/*							(ie/ every tree uses every edge and every via)				*/
/*		BAIL_ON_LAMBDA	:	We stop solving when our maximum congestion is below this	*/
/*							value.	(default = 1.0)										*/
/*		DELTA_LAM		:	We stop solving when the change in lamda is less than this.	*/
/*							(default = 0)												*/
/*		BENCH_STYLE		:	The type of benchmark for the input file					*/
/*							0 : IBM style (default)										*/
/*							1 : MCNC style												*/
/*							(these values are defined as macros)						*/
/*								IBM_BENCH, YAL_BENCH									*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/

#define HYBRID 0
#define CUTOFF 2


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "globals.h"

#include "graph.h"
#include "graphalg.h"
#include "list.h"
#include "set.h"
#include "helper.h"
#include "graphtype.h"
#include "linesearch.h"
#include "parseargs.h"
#include "rounding.h"
#include "routing.h"
#include "steiner.h"

#define LAYERS 2

#define QUERY_MODE 0
#define OUT_LEV 1

/* default parameter values */
#define N_DEFAULT 10
#define EPS_DEFAULT 0.001
#define TAU_HEURISTIC_DEFAULT 1
#define ALPHA_DEFAULT 0.5
#define BETA_DEFAULT 0.5
#define MULTIPLE_ROUND_DEFAULT 0
#define BAIL_ON_LAMBDA_DEFAULT 1.0
#define DELTA_LAM_DEFAULT 0
#define BENCH_STYLE_DEFAULT 0
#define CONGESTION_MAP_DEFAULT 0
#define WARM_START_DEFAULT 0

int main(int argc, char *argv[]) {

	/* graph data */
	Graph		grid,MLGraph;

	
	/* benchmark data*/
	int			width,height;		/* grid dimenstion*/
	int			vc,hc;				/* capacities */
	double		vl,hl;				/* physical length (given in MCNC benchmarks)*/
	Net			*T;					/* information about each net stored here*/
	CoordData	*unique_terminals;	/* set of unique terminals in the benchmark */
	int			no_unique_terminals;
	
	
	/* temporary counters*/
	int			temp_i,temp_j;
	double		partial;

	/* timer variables*/
	time_t		start,stop;
	double		ET;
	int			hours,mins,sec;
	time_t		iteration_start,iteration_stop;
	double		iteration_ET;
	
	/* sructure containing parameter settings*/
	args		the_parameters;
	
	
/* GLOBAL ROUTING VARIABLES----------------------------------------------------------------------- */

	double			EPS;		/* in (0,1], given accuracy*/
	int				N;			/* given bound on the number of iterations*/
	int				K;			/* number of nets*/
	double			*x;			/* (N+1)*K dimensional vector of indicator variables*/
	double			*f;			/* m dimensional vector, equivalent to edge congestion*/
	double			*f_bar;
	double			*p;			/* (m+1) dimensional price vector, in [0,1]*/
	double			SIG;		/* in (0,1], accuracy at current scaling phase*/
	double			LAM;		/* in R+, maximum component of f*/
	double			LAM_BAR;	
	double			THE;
	double			NU;			/* measure of duality gap*/
	double			TAU;		/* in (0,1], the step length*/
	double			*w;			/* an (m+n) dimensional vector, edge length of two layer graph H*/
	double			OMG;		/* scaling tolerence*/
	
	int				*c;			/* capacity for edge i (m+n) dimensional*/
	int				*c_orig;
	
	double			*L;			/* length of edge i in grid graph G (m dimensional)*/
	
	int				i;			/* index of edges (constraints)*/
	int				j;			/* index of trees (variables)*/
	int				k;			/* index of requests*/
	int				l;			/* index of iterations*/
	
	int				m;			/* Edges in grid graph*/
	int				n;			/* Vertices in grid graph*/
	
	double			g;			/* guessed objective value*/
	double			ALPHA,BETA;	/* artificial weights*/
	
	int				finished_scaling,finished_coordination;
	int				TAU_HEURISTIC;
	int				MULTIPLE_ROUND;
	double			BAIL_ON_LAMBDA;
	double			DELTA_LAM;
	int				BENCH_STYLE;
	double			old_lamda;
	int				CONGESTION_MAP;
	int				WARM_START;
	char			WS_FILE[32];
/* ------------------------------------------------------------------------------------------------*/	
	
	
	
	/* lookups*/
	AdjList****	Edge2Adj;
	
	/* debug*/
	double		avg_congestion;
	double		punishment;			

	double		total_cost,total_vias,objective_value;
			
	
	/* output variables*/
	int			retval;	
	FILE			*AllTrees;
	FILE			*Stats;
	int				num_routed;			/* never used*/
/*----------------------------------------------------------------------------------------------------------------------------------------
		END OF DECLARATIONS
----------------------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------
 SET PARAMETERS
------------------------------------------------------------------------------------------*/
	
	/* parse the argument file and set parameters*/
	parse_args(argv[2],&the_parameters);
	if (the_parameters.set_N)
		N = the_parameters.N;
	else
		N = N_DEFAULT;
	if (the_parameters.set_TAU)
		TAU = the_parameters.TAU;
	
	if (the_parameters.set_ALPHA)
		ALPHA = the_parameters.ALPHA;
	else
		ALPHA = ALPHA_DEFAULT;
	
	if (the_parameters.set_BETA)
		BETA = the_parameters.BETA;
	else
		BETA = BETA_DEFAULT;
		
	if (the_parameters.set_EPS)
		EPS = the_parameters.EPS;
	else
		EPS = EPS_DEFAULT;		
	
	if ( (the_parameters.set_TAU_HEURISTIC)&&(the_parameters.set_TAU == 0))
		TAU_HEURISTIC = the_parameters.TAU_HEURISTIC;
	else
		TAU_HEURISTIC = TAU_HEURISTIC_DEFAULT;

	if (the_parameters.set_MULTIPLE_ROUND)
		MULTIPLE_ROUND = 1;
	else
		MULTIPLE_ROUND = MULTIPLE_ROUND_DEFAULT;
	
	if (the_parameters.set_BAIL_ON_LAMBDA)
		BAIL_ON_LAMBDA = the_parameters.BAIL_ON_LAMBDA;
	else
		BAIL_ON_LAMBDA = BAIL_ON_LAMBDA_DEFAULT;

	if (the_parameters.set_DELTA_LAM)
		DELTA_LAM = the_parameters.DELTA_LAM;
	else
		DELTA_LAM = DELTA_LAM_DEFAULT;
	
	if (the_parameters.set_BENCH_STYLE)
		BENCH_STYLE = the_parameters.BENCH_STYLE;
	else
		BENCH_STYLE = BENCH_STYLE_DEFAULT;
	
	if (the_parameters.set_CONGESTION_MAP)
		CONGESTION_MAP = the_parameters.CONGESTION_MAP;
	else
		CONGESTION_MAP = CONGESTION_MAP_DEFAULT;
	
	if(the_parameters.set_WARM_START) {
		WARM_START = 1;
		strcpy(WS_FILE,the_parameters.WARM_START);
	}
	else
		WARM_START = WARM_START_DEFAULT;

	
	/* used to output congestion maps*/
	if (CONGESTION_MAP) system("mkdir output");

	srand( (unsigned)time(NULL) );
	
	if (argc < 3) {
		printf("Usage : route bench_file param_file\n");
		printf("\twhere bench file specifies the name of the benchmark file and param_file\n");
		printf("\tis the name of the file where the parameter settings are contained\n");
		exit(1);
	}

/*------------------------------------------------------------------------------------------
	PARSE BENCHMARK FILE
------------------------------------------------------------------------------------------*/
	if (BENCH_STYLE == IBM_BENCH) {	
		retval = parseBenchmark(argv[1], &K, &T, N, &width, &height,&hc,&vc,&num_routed);
		hl = 1.0;
		vl = 1.0;
	}
	else						/* MCNC style*/
		parse_yal_benchmark(argv[1],&K, &T, N, &width, &height, &hc, &vc, &hl, &vl);
		


	
	/* initialize unique terminals.   We won't always be using this so this might be unitialized otherwise */
	unique_terminals = NULL;
	no_unique_terminals = -1;
	find_unique_terminals(T,K, &unique_terminals, &no_unique_terminals); 

	if (QUERY_MODE) {
		query_mode(T,K);
		exit(0);
	}
	
	
	/* print benchmark statistics*/
	printf("\n%s\n",argv[1]);
	printf("\tWidth: %d, Height: %d\n",width,height);
	printf("\tNumber of Nets: %d\n",K);
	printf("\tVertical Capacity: %d\n",vc);
	printf("\tHorizontal Capacity: %d\n",hc);
	if (BENCH_STYLE == YAL_BENCH) {
		printf("\tVertical Length:   %f\n",vl);
		printf("\tHorizontal Length: %f\n",hl);
	}
	printf("\n");
	

	if ((Stats = fopen("Stats.txt","w")) == NULL) {
		printf("Couldn't open Stats.txt for writing\n");
		exit(1);
	}
	
	fprintf(Stats,"Statistics File\n");
	fprintf(Stats,"---------------\n\n");
	
	fprintf(Stats,"Benchmark file : %s\n",argv[1]);
	fprintf(Stats,"\tWidth %d , Height %d\n",width,height);
	fprintf(Stats,"\tVertical Capacity %d , Horizontal Capacity %d\n",vc,hc);
	if (BENCH_STYLE == YAL_BENCH) {
		fprintf(Stats,"\tVertical Length:   %f\n",vl);
		fprintf(Stats,"\tHorizontal Length: %f\n",hl);
	}
		
	fprintf(Stats,"\tNum Nets %d\n\n",K);
	



/*------------------------------------------------------------------------------------------
 CREATE GRAPHS
------------------------------------------------------------------------------------------*/
	gen_gridgraph(&grid,width,height);
	printf("Grid Graph: V = %d, E = %d\n",grid.vcount, grid.ecount);				/* create a 2d grid graph from the benchmark file*/
	fprintf(Stats,"Grid Graph: V = %d, E = %d\n",grid.vcount, grid.ecount);


	grid2VL(&grid,&MLGraph);														/* convert grid graph to multilayer graph (by adding */
	printf("2-Layer Graph: V = %d, E = %d\n",MLGraph.vcount, MLGraph.ecount);		/* vias and additional edges)*/
	fprintf(Stats,"2-Layer Graph: V = %d, E = %d\n",MLGraph.vcount , MLGraph.ecount);
		   
		   
	/*-----------------------------------------------------------
	 Create lookup arrays
	-----------------------------------------------------------*/	
	
	if ((Edge2Adj = (AdjList****) malloc(sizeof(AdjList***)*width))	== NULL) {	/* create function (Edge2Adj) from edges to adjacency lists*/
		printf("Edge2Adj : Memory allocation error\n"); fflush(stdout);
		exit(1);
	}
		
	for (temp_i = 0; temp_i < width; temp_i++)
		if ((Edge2Adj[temp_i] = (AdjList***) malloc(sizeof(AdjList**)*height))==NULL) {
			printf("Edge2Adj : Memory allocation error\n"); fflush(stdout);
			exit(1);
		}
	
	for (temp_i = 0; temp_i < width; temp_i++)
		for (temp_j = 0; temp_j < height; temp_j++)
			if ((Edge2Adj[temp_i][temp_j] = (AdjList**) malloc(sizeof(AdjList*)*LAYERS))==NULL){
				printf("Edge2Adj : Memory allocation error\n"); fflush(stdout);
				exit(1);
			}
	
	set_Edge2Adj(&MLGraph,Edge2Adj);	


	init_global_heap(width,height);
	init_shortest_path_globals(width,height,no_unique_terminals);
	init_xy2term(width,height,unique_terminals,no_unique_terminals);
/*
	for (i = 0; i < no_unique_terminals; i++) 
		printf("%d : (%d,%d)\n",i,unique_terminals[i].x,unique_terminals[i].y);
	printf("\n");

*/
		
/*---------------------------------------------------------------------------------------------------------------------------------------------
 INITIALIZATION
---------------------------------------------------------------------------------------------------------------------------------------------*/
		n = grid.vcount;
		m = grid.ecount / 2;					/* divide by two because each edge is two half edges*/
	
		if (the_parameters.set_g)				/* MUST HAVE VALUES FOR ALPHA,BETA,m,n*/
			g = the_parameters.g;
		else {
			double max;
			
			if (hl > vl)
				max = hl;
			else
				max = vl;
			
			g = K*(ALPHA*m*max + BETA*n);
		}
	
	
		
		punishment = 1.0;
		start = time(NULL);
	

		
		fprintf(Stats,"m = %d , n = %d\n",m,n);
		fprintf(Stats,"g = %f\n\n",g);
		fprintf(Stats,"TAU_HEURISTIC = %d\n",TAU_HEURISTIC);
		fprintf(Stats,"EPS = %E\n",EPS);
		fprintf(Stats,"N = %d\n",N);
		fprintf(Stats,"Alpha = %f, Beta = %f\n",ALPHA,BETA);
		
		fprintf(Stats,"Iteration\tSIG\t\tOMG\t\tTAU\t\tLAM\tLAM_BAR\t\tNU\t\tTHE\t\tOBJ\t\tTIME\n");
		fprintf(stdout,"Iteration\tSIG\t\tOMG\t\tTAU\t\tLAM\tLAM_BAR\t\tNU\t\tTHE\t\tOBJ\t\tTIME\n");
		
		
		
		if ((c = (int*) malloc(sizeof(int)*(m+1+1)))==NULL) {
			printf("c : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		set_capacity_lookup(&c,&MLGraph,hc,vc);
		c[m+1] = g;

		if (HYBRID) {
			if ((c_orig = (int*) malloc(sizeof(int)*(m+1+1)))==NULL) {
				printf("c_orig : memory allocation error\n"); fflush(stdout);
				exit(1);
			}
			set_capacity_lookup(&c_orig,&MLGraph,hc,vc);
			c_orig[m+1] = g;
		}

	
		
		if ((L = (double*) malloc(sizeof(double)*(m+1)))==NULL) {
			printf("L : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		
		if (BENCH_STYLE == YAL_BENCH)
			set_length_lookup(&L,&MLGraph,hl,vl);			
		else {
			for (i = 1; i <= m; i++) {
				L[i] = 1.0;
			}
		}
			
		
		if ((x = (double*) malloc(sizeof(double)*((N+1)*K + 1)))==NULL) {
			printf("x : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		
		
		for (j = 1; j <= K; j++)
			x[j] = 1.0;
		for (j = K+1; j <= (N+1)*K; j++)
			x[j] = 0.0;
				
		SIG = 1.0;
		
		OMG = (1.0 + SIG) / ( ( 1.0 + (SIG/6.0) )*(m+1.0) );
	
		if ((p = (double*) malloc(sizeof(double)*(m+2)))==NULL) {
			printf("p : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		for (i = 1; i <= m+1; i++) 
			p[i] = 1.0/(m+1.0);
		
		/* calculate weights */
		if ((w = (double*) malloc(sizeof(double)*(m+n+1)))==NULL) {
			printf("w : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		for (i = 1; i <= m; i++) {
			if (c[i] != 0.0)
				w[i] = punishment*(p[i]/c[i] + ALPHA*p[m+1]*L[i]/g);
			else
				w[i] = DBL_MAX;
		}
		for (j = 1; j <= n; j++)
			w[m+j] = p[m+1]*BETA/g;

				
		sync_weights(&MLGraph,w);	
		
															/* sets all edges in MLGraph to weights defined in w array*/
		iteration_start = time(NULL);


		/* generate trees for nets with degree less than or equal to CUTOFF */
		if (HYBRID){										
			
			for (k = 1; k <= K; k++) {
				
				if (T[k].no_terminals <= CUTOFF) {
					
					gridSteinerFH(&MLGraph, Edge2Adj, width, height, T[k].terminals, T[k].no_terminals, L, T[k].net_num, 
								  &(T[k].edge_count), &(T[k].SteinerTree),0,K);
					
					for (l = 0; l <= N; l++)
						fix_net(T,l*K + k);

					T[k].tree_cost = 0.0;
					for (i = 0; i < T[k].edge_count; i++)
						if (T[k].SteinerTree[i] <= m) 
							T[k].tree_cost += L [ T[k].SteinerTree[i] ];
						
					
					T[k].vias = 0;
					for (i = 0; i < T[k].edge_count; i++)
						if (T[k].SteinerTree[i] > m)
							T[k].vias++;
					
					T[k].routed = 1;
				}
			}

			/* adjust capacities */
			for (k = 1; k <= K; k++) {
				if (T[k].no_terminals <= CUTOFF) {
					int remove;
					
					remove = 0;
					for (i = 0; i < T[k].edge_count; i++) {
						if (T[k].SteinerTree[i] <= m)
							c[ T[k].SteinerTree[i] ]--;
						
						/* bail!! */
						if (c[ T[k].SteinerTree[i] ] == 0) 
							remove = 1;
					}
					
					if (remove) {
						
						/* re-adjust capacity back to the way it was */
						for (i = 0; i < T[k].edge_count; i++) 
							if (T[k].SteinerTree[i] <= m)
								c[ T[k].SteinerTree[i] ]++;
						
						T[k].fixed = 0;
					
					}
				}
			}
		}


	
		if (OUT_LEV>1) printf("Generating first set of trees\n");
		if (!WARM_START) {
			for (k = 1; k <= K; k++) {
				
				if  ( !T[k].fixed  || (!HYBRID) ) {
					gridSteinerFH(&MLGraph, Edge2Adj, width, height, T[k].terminals, T[k].no_terminals, L, T[k].net_num, 
								  &(T[k].edge_count), &(T[k].SteinerTree),0,K);
							
		/*			if (T[k].no_terminals <= 3)
					gen_tree_FLUTE(&MLGraph, Edge2Adj, T[k].terminals, T[k].no_terminals, &(T[k].SteinerTree), &(T[k].tree_cost));
		*/						
					/* compute cost*/
					T[k].tree_cost = 0.0;
					for (i = 0; i < T[k].edge_count; i++)
						if (T[k].SteinerTree[i] <= m)
							T[k].tree_cost += L [ T[k].SteinerTree[i] ];
					
					T[k].vias = 0;
					for (i = 0; i < T[k].edge_count; i++)
						if (T[k].SteinerTree[i] > m)
							T[k].vias++;
							
					T[k].routed = 1;
				}
				
				if (OUT_LEV > 1) {			  
					if ( k == K/4 )
						printf("\t25%d\n",(char)37);
					if (k == K/2 )
						printf("\t50%d\n",(char)37);
					if (k == (K/4 + K/2))
						printf("\t75%d\n",(char)37);
				}
			}
		}
		else {
			FILE *trees;
			
			trees = fopen(WS_FILE,"r");
			fscanf(trees,"%*s");				/* #######################  */
				   
			   for (k = 1; k <= K; k++)
			   {
				   int		edge_count;
				   int		i,index,found;
				   double	tree_cost;
				   int	*SteinerTree;
				   char		ID[50];
				   
				   get_next_tree(trees,&MLGraph,hl,vl,&SteinerTree,&edge_count,&tree_cost,ID);
				   
				   for (found = 0, i = 1; ((i<=K)&&(!found)); i++)
					   if ( !strcmp(ID,T[i].ID) ) {
						   found = 1;
						   index = i;
					   }
					
				   if (!found) {
					   printf("ERROR: not found.   k = %d, ID=%s\n",k,ID);
					   exit(1);
				   }
				   
				   
				   T[index].edge_count = edge_count;
				   T[index].SteinerTree = SteinerTree;
				   T[index].routed = 1;
				   T[index].tree_cost = tree_cost;
				   T[index].vias = 0;
				   T[index].fixed = 0;
				}
		}


//printf("\n\n");
//for (k = 1; k <= K; k++) {
//	int i,count;
//	if (T[k].fixed)
//		printf("*");
//	else
//		printf(" ");
//	for (count = 0,i = k; count <= l; i+=K,count++) 
//		printf("x[%2d] = %f\t",i,x[i]);
//	printf("\n");
//}
//printf("\n\n");

		reset_shortest_path(no_unique_terminals);
		

		if ((f = (double*) malloc(sizeof(double)*(m+1+1)))==NULL) {
			printf("f : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		
		if ((f_bar = (double*) malloc(sizeof(double)*(m+1+1)))==NULL) {
			printf("f_bar : memory allocation error\n"); fflush(stdout);
			exit(1);
		}
		
		avg_congestion = 0.0;
		for (i = 1; i <= m; i++) {

			f[i] = 0.0;			
			for (j = 1; j <= K; j++)
				f[i] += (double)U(T,i,j);
			
			f[i] = f[i]/(double)c[i];
			avg_congestion += f[i];		
		}
		avg_congestion /= m;
		
		f[m+1] = 0.0;
		for (j = 1; j <= K; j++)
			f[m+1] += T[j].tree_cost;
		f[m+1] *= ALPHA/g;
		
		if (CONGESTION_MAP) {
			edge_Histogram(&MLGraph,f_bar,1,m,"Congestion_iter",0);
			edge_Histogram(&MLGraph,f,1,m,"Congestion",0);
			system("mv Congestion* output");
		}
		
		partial = 0.0;
		for (j = 1; j <= K; j++)
			partial += (double)T[j].vias;
		partial *= BETA/g;
		
		f[m+1] += partial;
		
		if (OUT_LEV > 1) printf("Average congestion after initialization = %E\n",avg_congestion);
				
		LAM = -DBL_MAX;
		for (i = 1; i <= m+1; i++)
			if (f[i] > LAM)
				LAM = f[i];
		
		LAM_BAR = LAM;
		old_lamda = LAM;

		
		iteration_stop = time(NULL);
		iteration_ET = difftime(iteration_stop,iteration_start);
		hours = ((int)iteration_ET) / 3600;
		mins  = ((int)iteration_ET - (3600*hours)) / 60;
		sec   = (int)iteration_ET - (3600*hours) - 60*mins;
		
		
		fprintf(Stats,"0\t%E\t%E\t\t\t%E\t%E\t\t\t\t\t\t\t%d,%d,%d\n",SIG,OMG,LAM,LAM_BAR,hours,mins,sec);
		fprintf(stdout,"0\t%E\t%E\t\t\t%E\t%E\t\t\t\t\t\t\t%d,%d,%d\n",SIG,OMG,LAM,LAM_BAR,hours,mins,sec);
			   
		finished_scaling = 0;

		l = 0;

		if (OUT_LEV>1)printf("Start scaling...\n");
		while ((!finished_scaling)&&(l<N)) {
			
			finished_coordination = 0;
			
			if (OUT_LEV>1)printf("\tStart coordination...\n");
			while ((!finished_coordination)&&(l < N)) {
				double	LB,UB,Tol;
				int		integer_feaseable;
				
				l++;
				if (OUT_LEV>1) printf("\tEntering iteration %d\n",l);
				iteration_start = time(NULL);

				LB = LAM / (1.0 - (SIG/6.0)*(1.0/(m+1.0)));
				UB = LAM / (1.0 - (SIG/6.0));
				Tol = (EPS*EPS)/(m+1.0);
				
				
				if (OUT_LEV>1) printf("\t\tSolve for theta..."); fflush(stdout);
				THE = solve_BISECT_THETA(LB, UB , Tol, SIG,m,f,MAX_BISECT_ITERS);
				if (OUT_LEV>1) printf("done\n"); fflush(stdout);
				
				
				for (i = 1; i <= m+1; i++) {
					p[i] = (SIG*THE) / ( 6.0*(m+1)*(THE-f[i]) );
				}
				
				for (i = 1; i <= m; i++) {
					if (c[i] != 0.0)
						w[i] = punishment*(p[i]/c[i] + ALPHA*p[m+1]*L[i]/g);
					else
						w[i] = DBL_MAX;
				}
				
				for (j = 1; i <= n; i++) {
					w[m+j] = p[m+1]*BETA/g;
				}
				
				sync_weights(&MLGraph,w);
								
				if (OUT_LEV>1) printf("\t\tGenerating trees for iteration %d...\n",l);
				for (k = 1; k <= K; k++) {
					
					if (!T[k].fixed) {
					
						gridSteinerFH(&MLGraph, Edge2Adj, width, height, T[l*K + k].terminals, T[l*K + k].no_terminals, 
						L, T[k].net_num, &(T[l*K + k].edge_count), &(T[l*K + k].SteinerTree),l,K);
											
						/* compute cost*/
						T[l*K + k].tree_cost = 0.0;
						for (i = 0; i < T[l*K + k].edge_count; i++)
							if (T[l*K + k].SteinerTree[i] <= m)
								T[l*K + k].tree_cost += L [ T[l*K + k].SteinerTree[i] ];
						
						T[l*K + k].vias = 0;
						for (i = 0; i < T[l*K + k].edge_count; i++)
							if (T[l*K + k].SteinerTree[i] > m)
								T[l*K + k].vias++;
						
						T[l*K + k].routed = 1;
						
						if (OUT_LEV>1) {
							if ( k == K/4 )
								printf("\t25%d\n",(char)37);
							if (k == K/2 )
								printf("\t50%d\n",(char)37);
							if (k == (K/4 + K/2))
								printf("\t75%d\n",(char)37);
						}
					}
				}
				
				reset_shortest_path(no_unique_terminals);


				if (OUT_LEV>1) printf("\t\tdone.\n");				
				fflush(stdout);
				
				
				for (i = 1; i <= m; i++) {
					
					partial = 0.0;
					for (j = l*K + 1; j <= (l+1)*K; j++)
						partial += (double) U(T,i,j);
					
					if (c[i] != 0.0)
						f_bar[i] = partial/c[i];
					else
						f_bar[i] = DBL_MAX;
				}

				
				partial = 0.0;
				for (j = l*K + 1; j <= (l+1)*K; j++)
					partial += T[j].tree_cost;
				
				f_bar[m+1] = partial*ALPHA/g;
				
				partial = 0.0;
				for (j = l*K + 1; j <= (l+1)*K; j++)
					partial += (double)T[j].vias;
			
				f_bar[m+1] += partial*BETA/g;
				
				/*--------------------------------------------------------------------------*/
				/*  Checking integer feaseability											*/
				/*--------------------------------------------------------------------------*/
				integer_feaseable = 1;
				for (i = 1; i <= m+1; i++) 
					if (f_bar[i] >= 1.0) integer_feaseable = 0;
				if (integer_feaseable)
					if (OUT_LEV>1) printf("********** ITERATION %d :  INTEGER FEASEABLE **********\n",l);
				
				partial = 0.0;
				for (i = 1; i <= m+1; i++)
					partial += p[i]*(f[i] - f_bar[i]);
				NU = partial;
				
				
				partial = 0.0;
				for (i = 1; i <= m+1; i++)
					partial += p[i]*(f[i] + f_bar[i]);

				NU /= partial;
				
/*				if ( NU <= (SIG/6.0) ) {*/
				if ( 0 ) {
					finished_coordination = 1;
					if (OUT_LEV>1) printf("\tFinished coordination\n\n**************************");
					if (OUT_LEV>1) printf("\t\tNU = %f , SIG/6.0 = %f\n",NU,SIG/6.0);
				}
				else {
					
					TAU = SIG*THE*NU / (12.0*(m+1)*partial);			/* DEFAULT*/

					if (TAU_HEURISTIC == 1) 
						TAU = SIG*THE*NU / (12.0*(m+1)*partial);		/* HEURISTIC 1*/
					
					if (TAU_HEURISTIC == 2) 
						TAU = 1.0 - NU;									/* HEURISTIC 2*/

					if (TAU_HEURISTIC == 3) 
						TAU = NU;										/* HEURISTIC 3*/

					if (TAU_HEURISTIC == 5) 
						TAU = ((1+sqrt(5))/2)-1;						/* HEURISTIC 5*/
					
					if (TAU_HEURISTIC == 6) 
						TAU = 0.040;									/* HEURISTIC 6*/
					
					if (TAU_HEURISTIC == 7) {							/* HEURISTIC 7*/
						if (NU > 1-NU)
							TAU = NU;
						else
							TAU = 1 - NU;
					}
					
					if (TAU_HEURISTIC == 8) {					/* HEURISTIC 8 : Line Search to minimize violation*/
						double tmp_f[m+2];
						double interval;
						double bestTau,minval;
						double violation;
					
						minval = DBL_MAX;								
						bestTau = 0.0;
						for (interval = 0.0001; interval < 0.5; interval = interval + 0.0001) { 
							
							/* create a copy of f[i]*/
							for (i = 1; i <= m+1; i++)
								tmp_f[i] = f[i];

							/* calculate updated LHS*/
							for (i = 1; i <= m+1; i++)
								tmp_f[i] = (1.0 - interval)*f[i] + interval*f_bar[i];
							
							/* calculate violation*/
							violation = 0.0;
							for (i = 1; i <= m+1; i++)
								if (tmp_f[i] > 1.0)
									violation += tmp_f[i] - 1.0;
							
							if (violation < minval) {
								minval = violation;
								bestTau = interval;
							}
						}
						TAU = bestTau;
					}
					
					if (TAU_HEURISTIC == 9) {					/*HEURISTIC 9 : Line search to minimize potentiol function*/
						double *tmp_f;
						double interval;
						double bestTau,minval;
						double pot_val;
						int		found_neg;
						
						minval = DBL_MAX;
						bestTau = 0.0;
						if ((tmp_f = (double*) malloc(sizeof(double)*(m+2)))==NULL) {
							printf("main.c : tmp_f mem allocation error\n"); fflush(stdout);
							exit(1);
						}
						
						
						minval = DBL_MAX;
						found_neg = 0;
						for (interval = 0.0001; interval < 0.4; interval = interval + 0.0001) { 
							double theta,tmp_LAM;
							
							for (i = 1; i <= m+1; i++)
								tmp_f[i] = f[i];
							
							
														
							/* calculate updated LHS*/
							for (i = 1; i <= m+1; i++)
								tmp_f[i] = (1.0 - interval)*f[i] + interval*f_bar[i];
							
							tmp_LAM = -DBL_MAX;
							for (i = 1; i <= m+1; i++)
								if (tmp_f[i] > tmp_LAM)
									tmp_LAM = tmp_f[i];
							
							
							
							/* solve for theta*/
							LB = tmp_LAM / (1.0 - (SIG/6.0)*(1.0/(m+1.0)));
							UB = tmp_LAM / (1.0 - (SIG/6.0));
							Tol = (EPS*EPS)/(m+1.0);							
							theta = solve_BISECT_THETA(LB, UB , Tol, SIG,m,tmp_f,MAX_BISECT_ITERS);
							
							/* calculate potential value*/
							pot_val = log(theta);
							partial = 0.0;
							for (i = 1; i <= m+1; i++)
								partial += log(theta - tmp_f[i]);
							partial *= ((SIG/6.0)/(m+1.0));
							pot_val = pot_val - partial;
							
							
							if (pot_val < minval) {
								minval = pot_val;
								bestTau = interval;
							}
						
						}
						TAU = bestTau;
						free(tmp_f);
					}
					
					if (TAU_HEURISTIC == 10) {			/* secant method to minimize potential function*/
						double x_0,x_1;
						double epsilon,delta;
						int	   max_iters;
						double *tmp_f;
						
						epsilon = 0.001;
						delta =	epsilon;
						max_iters = 50;
						x_0 = 0.2;
						x_1 = 0.19;
						
						/* create a workspace*/
						if ((tmp_f = (double*) malloc(sizeof(double)*(m+2)))==NULL) {
							printf("main.c : tmp_f mem allocation error\n"); fflush(stdout);
							exit(1);
						}
						
						TAU = secant(x_0,x_1,epsilon,delta,max_iters,pot_fun,f,f_bar,tmp_f,SIG,EPS,m);
						
						free(tmp_f);
					}
					
					if (TAU_HEURISTIC == 11) {
						double a,b,delta;
						double *tmp_f;
						
												
						if ((tmp_f = (double*) malloc(sizeof(double)*(m+2)))==NULL) {
							printf("main.c : tmp_f mem allocation error\n"); fflush(stdout);
							exit(1);
						}
						a = 0.0;
						b = 1.0;
						delta = 0.00001;
						
						TAU = solve_BISECT_potfun(a,b,delta,delta,SIG,m,f_bar,f,tmp_f,EPS,MAX_BISECT_ITERS);
						
						if (TAU == 0.0) {
							TAU_HEURISTIC = 1;
							TAU = SIG*THE*NU / (12.0*(m+1)*partial);		/* HEURISTIC 1*/
						}
						free(tmp_f);
					}
							
										
							
					for (j = 1; j <= l*K; j++) {
						if (T[j].fixed)
							continue;
						else 
							x[j] = (1-TAU)*x[j];
					}
					
					
					for (j = l*K + 1; j <= (l+1)*K; j++) {
						if (T[j].fixed)
							x[j] = 0.0;
						else
							x[j] = TAU;
					}

//printf("\n\n");
//for (k = 1; k <= K; k++) {
//	int i,count;
//	if (T[k].fixed)
//		printf("*");
//	else
//		printf(" ");
//	for (count = 0,i = k; count <= l; i+=K,count++) 
//		printf("x[%2d] = %f\t",i,x[i]);
//	printf("\n");
//}
//printf("\n\n");

					avg_congestion = 0.0;
					for (i = 1; i <= m+1; i++) {
						f[i] = (1-TAU)*f[i] + TAU*f_bar[i];
						avg_congestion += f[i];
					}
					avg_congestion /= m+1;
					
					if (CONGESTION_MAP) {
						edge_Histogram(&MLGraph,f_bar,1,m,"Congestion_iter",l);
						edge_Histogram(&MLGraph,f,1,m,"Congestion",l);
						system("mv Congestion* output");
					}
						
					
					LAM = -DBL_MAX;
					for (i = 1; i <= m+1; i++)
						if (f[i] > LAM)
							LAM = f[i];
				}
				
				if ( (LAM <= OMG*LAM_BAR)||(NU < SIG/6.0) ) {
					finished_coordination = 1;
					if (OUT_LEV > 1) printf("\tFinished coordination\n\n**************************");
					if (OUT_LEV > 1) printf("\tLAM = %E less than OMG*LAM_BAR = %E\n",LAM,OMG*LAM_BAR);
				}
				
				if (( LAM <= BAIL_ON_LAMBDA )||(NU <= 0)) {			/* quit if max congested edge is below this value or if NU has gone negative*/
					finished_coordination = 1;
					finished_scaling = 1;
				}
				
				if ( fabs(LAM - old_lamda) < DELTA_LAM ) {
					finished_coordination = 1;
					finished_scaling = 1;
				}
				else 
					old_lamda = LAM;
				
				total_vias = 0.0;
				total_cost = 0.0;
				for (k = 1; k <= (l+1)*K; k++) {
					total_vias += x[k]*T[k].vias;
					total_cost += x[k]*T[k].tree_cost;
				}
				objective_value = ALPHA*total_cost + BETA*total_vias;
				
				iteration_stop = time(NULL);
				iteration_ET = difftime(iteration_stop,iteration_start);
				hours = ((int)iteration_ET) / 3600;
				mins  = ((int)iteration_ET - (3600*hours)) / 60;
				sec   = (int)iteration_ET - (3600*hours) - 60*mins;
				
					
				fprintf(Stats,"%d\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%d,%d,%d\n",l,SIG,OMG,TAU,LAM,LAM_BAR,NU,THE,objective_value,hours,mins,sec);
				fprintf(stdout,"%d\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%d,%d,%d\n",l,SIG,OMG,TAU,LAM,LAM_BAR,NU,THE,objective_value,hours,mins,sec);

			}	
			
			SIG = SIG/2.0;
			LAM_BAR = LAM;
			OMG = (1+SIG)/(1+(2.0*SIG));
			
			if (SIG <= EPS/2) {
				finished_scaling = 1;
				if (OUT_LEV>1) printf("Finished scaling\n");
			}
		}
			
		if (OUT_LEV>1) printf("DONE!\n");
		

		/* Print the Fractional Solution*/
		if ((AllTrees = fopen("AllTrees.txt","w"))!=NULL) { 
			print_all_nets(AllTrees,T,N,K);
			fclose(AllTrees);
		}
		else 
			printf("Couldn't open AllTrees.txt for writing\n");			
		

		/* print the statistics file*/
		if (Stats != NULL) {
			double	objective_value;
			double	total_cost;
			double	total_vias;			
			int		feaseable;
			int		infeaseable_edges[m+2];
			
			
			feaseable = 1;
			total_vias = 0.0;
			total_cost = 0.0;
			for (k = 1; k <= (l+1)*K; k++) {
				total_vias += x[k]*T[k].vias;
				total_cost += x[k]*T[k].tree_cost;
			}
			objective_value = ALPHA*total_cost + BETA*total_vias;

			if (!HYBRID) {

				fprintf(Stats,"Checking feaseability of the (fractional) solution:\n");
				for (i = 1; i <= m+1; i++) {
					infeaseable_edges[i] = 0;
					if (f[i] > 1.0) {
						infeaseable_edges[i] = 1;
						feaseable = 0;
					}
				}
						
				if (feaseable)
					fprintf(Stats,"\tFeaseable!\n");
				else {
					fprintf(Stats,"\tInfeaseable\n");
					fprintf(Stats,"\tConstraint Violations:\n");
					for (i = 1; i <= m+1; i++) 
						if (infeaseable_edges[i]) 
							fprintf(Stats,"\t\tEdge %d: Used %f times, capacity is %d\n",i,f[i]*c[i],c[i]);
				}
			}
			
			fprintf(Stats,"Obj Val = %e\n",objective_value);
			fprintf(Stats,"\tAlpha = %f , Beta = %f, total tree cost = %f , total_vias = %f\n",ALPHA,BETA,total_cost,total_vias);

			printf("Objective Value of FRACTIONAL solution = %e\n",objective_value);
			printf("\tAlpha = %f , Beta = %f, total tree cost = %f , total_vias = %f\n\n",ALPHA,BETA,total_cost,total_vias);
							
			
		}
				
		
			
		
			
/*------------------------------------------------------------------------------------------
 CLEAN-UP
------------------------------------------------------------------------------------------*/
    /* calculate time taken */
	stop = time(NULL);
	ET = difftime(stop,start);
	hours = ((int)ET) / 3600;
	mins  = ((int)ET - (3600*hours)) / 60;
	sec   = (int)ET - (3600*hours) - 60*mins;
	printf("TOTAL TIME for file %s ---  %d:%d:%d\n",argv[1],hours,mins,sec);
	fprintf(Stats,"TOTAL TIME for file %s ---  %d:%d:%d\n",argv[1],hours,mins,sec);
	fclose(Stats);

	for (i = 0; i < width; i++)									/* free lookups*/
		for (j = 0; j < height; j++) 
			free(Edge2Adj[i][j]);
	for (i = 0; i < width; i++)
		free(Edge2Adj[i]);	
	free(Edge2Adj);
	
	/* free terminal sets, for each iteration, the terminals point to the same memory */
	for (i = 1; i <= K; i++) {								/* free output*/
		free(T[i].terminals);
	}
	
	/* need to free each one because all trees have their own memory, unlike terminal sets */
	for (i = 1; i <= K*l; i++) 
		free(T[i].SteinerTree);
	free(T);
	
	free(c);
	free(L);
	free(x);
	free(p);
	free(w);
	free(f);
	free(f_bar);
	
	if (unique_terminals != NULL)
		free(unique_terminals);
	destroy_global_heap();
	destroy_shortest_path_globals(width,height,no_unique_terminals);
	graph_destroy(&grid);										/* deystroy the graphs*/
    graph_destroy(&MLGraph);

	return 0;
}
