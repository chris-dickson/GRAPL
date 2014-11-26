#define MAX_TRIES 10000

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "list.c"
#include "set.c"
#include "graph.c"
#include "graphalg.h"
#include "rounding.c"


#define MAX_WORD_SIZE 100
#define LINE_SIZE 1000	
#define MULTIPLE_ROUND 1

typedef struct CoordData_ {
	int x,y,z;
}CoordData;

typedef struct Net_ {
	char		num[10];
	int			net_num;
	int			no_terminals;
	CoordData	*terminals;
	int			*SteinerTree;
	int			edge_count;
	int			vias;
	double		tree_cost;
	int			routed;		// boolean, 0 if not routed
}Net;

#include "aux.c"

int main(int argc, char *argv[]) {
	FILE	*inSolutionFile, *inStatFile, *outFile;
	char	a_word[MAX_WORD_SIZE],a_line[LINE_SIZE];
	int		i,done,iters,l,k,j;
	
	int		Mrd;
	int		best_Mrd;
	
	int		N,K;
	double	*x;
	int		*x_star;
	double	*f;
	double	TAU;
	Net		*T,*Solution;
	
	Graph	G,H;
	int		width,height;
	int		hc,vc;
	double	vl,hl;
	int		*c;
	double	*L;
	int		m;
	
	double	ALPHA,BETA;
	char	float_value[32];
	double	total_cost,total_vias;
	double	objective;
	double	best_objective;
	
	time_t		iteration_start,iteration_stop;
	double		iteration_ET;
	int			hours,mins,sec;
	
	
	if (argc != 4) {
		printf("Usage:  fractional_rounder solutionfile statfile output\n");
		printf("\t solutionfile : the fractional solution output file from route\n");
		printf("\t statfile : the statistics output file from route\n");
		printf("\t output : the solution file to be output\n");
		exit(1);
	}	
	
	//open input file
	inSolutionFile = fopen(argv[1],"r");
	inStatFile = fopen(argv[2],"r");
//	outFile = fopen(argv[3],"w");

	printf("How many iterations?\n");
	scanf("%d",&iters);
	
	printf("Width: ");
	scanf("%d",&width);
	
	printf("Height: ");
	scanf("%d",&height);
	
	printf("Vertical Capacity: ");
	scanf("%d",&vc);
	
	printf("Horizontal Capacity: ");
	scanf("%d",&hc);
	
//	printf("Vertical Length: ");
//	scanf("%f",&vl);
//	
//	printf("Horizontal Length: ");
//	scanf("%f",&hl);
	
	printf("Alpha: ");
	scanf("%s",float_value);
	ALPHA = atof(float_value);
	
	printf("Beta: ");
	scanf("%s",float_value);
	BETA = atof(float_value);
	
	
	printf("ECHO\n");
	printf("\tAlpha = %E\n",ALPHA);
	printf("\tBeta = %E\n",BETA);
	
	srand( time(NULL) );
	

	// figure out K
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Num")) {
			done = 1;
		}
	}
	fscanf(inStatFile,"%*s %d",&K);
	
	// figure out N
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"N")) {
			fscanf(inStatFile,"%s",a_word);
			if (!strcmp(a_word,"=")) {
				done = 1;
			}
		}
	}
	fscanf(inStatFile,"%d",&N);
	
	x = (double*) malloc(sizeof(double)*((N+1)*K));


	if ((x_star = (int*) malloc(sizeof(int)*((N+1)*K + 1)))==NULL) {
		printf("x_star : memory allocation error\n"); fflush(stdout);
		exit(1);
	}
	

	// find first line
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Iteration")) {
			done = 1;
		}
	}
	fgets(a_line,LINE_SIZE,inStatFile); // rest of first line
	fgets(a_line,LINE_SIZE,inStatFile); // second line, TAU = 1.0
	
	for (i = 1; i <= K; i++)
		x[i] = 1.0;

	done = 0;
	while (!done) {
		char *ret;
		
		ret = fgets(a_line,LINE_SIZE,inStatFile);
		if (ret == NULL)							// END OF FILE
			done = 1;
		else {
			int iteration;
			char TAU[MAX_WORD_SIZE];
			double TAU_VAL;
			int		k,j;
			
//			printf("%s",a_line);					// echo
//			fflush(stdout);
			sscanf(a_line,"%d %*s %*s %s %*s %*s %*s %*s %*s %*s\n",&iteration,TAU);
			TAU_VAL = atof(TAU);
//			printf("\tIteration = %d\n",iteration);
//			printf("\tTAU = %E\n\n",TAU_VAL);

			for (l = 0; l < iteration; l++)
				for (k = 1; k <= K; k++)
					x[l*K + k] = x[l*K + k]*(1-TAU_VAL);
			
			for (k = 1; k <= K; k++)
				x[iteration*K + k] = TAU_VAL;
			
			if (iteration == iters)
				done = 1;
		}
	}
	
	// We now have our x array.   Now we must parse the solution file before we attempt to round
	if ((T = (Net*) malloc(sizeof(Net)*((iters+1)*K + 1)))==NULL) {
		printf("benchmark, temp : memory allocation error\n"); fflush(stdout);
		exit(1);
	}

		
	
	// parse the trees
	for (l = 0; l <= iters; l++) {
		int iter;
		
		fscanf(inSolutionFile,"%*s %d",&iter);		// Iteration %d
		fscanf(inSolutionFile,"%*s");				// -------------------
		fscanf(inSolutionFile,"%*s");				// \n
		
		// for each tree
		for (k = 1; k <= K; k++) {
			char no_term_string[MAX_WORD_SIZE], coord[MAX_WORD_SIZE];
			
			fscanf(inSolutionFile,"%*s %*d");						// Net %d
			fscanf(inSolutionFile,"%s",no_term_string);			// no_terminals=%d
			i = 0;									// skip to the = sign
			while (no_term_string[i] != '=')
				i++;
			i++;									// point to the number
			T[l*K+k].no_terminals = atoi( &(no_term_string[i]) );	// convert
	
			T[l*K+k].terminals = (CoordData*) malloc(sizeof(CoordData)*T[l*K+k].no_terminals);
	

			fscanf(inSolutionFile,"%*s %*s");						// terminal set
			for (i = 0; i < T[l*K+k].no_terminals; i++) {
				int j;
				
				fscanf(inSolutionFile,"%*d %*s %s",coord);		// %d : (%d,%d,%d)
				sscanf(coord,"(%d,%d,%d)",&(T[l*K+k].terminals[i].x),&(T[l*K+k].terminals[i].y),&(T[l*K+k].terminals[i].z));				
			}
			
			// get edge count
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].edge_count = atoi(&a_line[i]);
			
			// get vias
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].vias = atoi(&a_line[i]);
			
			// get tree cost
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].tree_cost = atof(&a_line[i]);
			
			fgets(a_line,LINE_SIZE,inSolutionFile);		//junk
			fgets(a_line,LINE_SIZE,inSolutionFile);		// steiner tree
			
			T[l*K+k].SteinerTree = (int*) malloc(sizeof(int)*T[l*K+k].edge_count);
			

			for (i = 0; i < T[l*K+k].edge_count; i++) 
				fscanf(inSolutionFile,"%d",&(T[l*K+k].SteinerTree[i]));
			
			qsort(T[l*K+k].SteinerTree , T[l*K+k].edge_count , sizeof(int), &qsort_match);
				

			
			fgets(a_line,LINE_SIZE,inSolutionFile);			
			fgets(a_line,LINE_SIZE,inSolutionFile);			
			
			
		}
	}
	
			
	// We now have our fractional solution stored in T as it would be in route.   Now, we must make the graphs
	fclose(inSolutionFile);
	fclose(inStatFile);
	
	gen_gridgraph(&G,width,height);
	grid2VL(&G,&H,width,height);
	if ((c = (int*) malloc(sizeof(int)*(((G.ecount)/2)+1+1)))==NULL) {
		printf("c : memory allocation error\n"); fflush(stdout);
		exit(1);
	}
	set_capacity_lookup(&c,&H,hc,vc);
	
//	if ((L = (double*) malloc(sizeof(int)*(((G.ecount)/2)+1+1)))==NULL) {
//		printf("L : memory allocation error\n"); fflush(stdout);
//		exit(1);
//	}
//	set_length_lookup(&L,&H,hl,vl);

	printf("\nPerforming randomized rounding\n");	
	if ((Solution = (Net*) malloc(sizeof(Net)*(K+1)))==NULL) {
		printf("main.c : Solution mem allocation error\n");
		exit(1);
	}
	
	m = G.ecount / 2;
	
	l = iters;
	
	if (MULTIPLE_ROUND) {
		int edge_count[m+2];
		int violations;
		int best_violations;
		int	tries;
		int	x_star_inverse[K+1];		
//		double	total_cost,total_vias;
		
		best_objective = DBL_MAX;
		best_violations = m;
		violations = m;
		best_Mrd = m;
		
		tries = 1;
		iteration_start = time(NULL);
		while (tries < MAX_TRIES) {
			
			// perform a randomized rounding
			randomized_round(x,K,l,x_star);
			
			for (k = 1; k <= K; k++){
				int index_of_one,found;
				int j;
				
				for (j = 0, found = 0; (j <= l)&&(!found); j++) {
					if (x_star[j*K+k] == 1) {
						found = 1;
						index_of_one = j*K+k;
						x_star_inverse[k] = index_of_one;
					}

				}
				
			}
			
			violations = 0;
			for (i = 1; i <= m; i++)
				edge_count[i] = 0;
			
			// compute congestion
			for (k = 1; k <= K; k++) 
				for (i = 1; i <= m; i++) 
					edge_count[i] += U_sorted(T,i,x_star_inverse[k]);
			
			// check for edge violations
			Mrd = 0;
			for (i = 1; i <= m; i++) {
				if ((double)edge_count[i] > c[i])
					violations++;
				if (edge_count[i] > Mrd)
					Mrd = edge_count[i];
			}
			
			// compute objective value
			total_cost = 0.0;
			total_vias = 0.0;
			for (k = 1; k <= K; k++) {
				total_cost += (double) T[x_star_inverse[k]].tree_cost;
				total_vias += (double) T[x_star_inverse[k]].vias;
			}
			objective = ALPHA*total_cost + BETA*total_vias;			
			
			
			// is this round any better?
			if ((Mrd < best_Mrd)||((Mrd == best_Mrd)&&(objective < best_objective))) {

				// found a solution, stop the timer
				iteration_stop = time(NULL);
				iteration_ET = difftime(iteration_stop,iteration_start);
				hours = ((int)iteration_ET) / 3600;
				mins  = ((int)iteration_ET - (3600*hours)) / 60;
				sec   = (int)iteration_ET - (3600*hours) - 60*mins;
				
				if (Mrd < best_Mrd) {
					best_Mrd = Mrd;
					best_objective = DBL_MAX;
					printf("Number of violations = %d, Mrd = %d  (After %d tries and %dh %dm %ds)\n",violations,Mrd,tries,hours,mins,sec);
				}
				
				if (objective < best_objective) {
					best_objective = objective;
					printf("\tObjective = %E (After %d tries and %dh %dm %ds)\n",objective,tries,hours,mins,sec);
				

					if ((outFile = fopen(argv[3],"w"))!=NULL) {
						
						// overwrite old solution
						for (k = 1; k <= K; k++)
							memcpy( &Solution[k] , &T[x_star_inverse[k]] , sizeof(T[x_star_inverse[k]]) );
						
						// print header
						if (violations > 0) fprintf(outFile,"INFEASEABLE INTEGER SOLUTION\n");
						else fprintf(outFile,"INTEGER SOLUTION\n");
						fprintf(outFile,"----------------------------\n");
						
						fprintf(outFile,"OBJECTIVE VALUE = %E\n\tWirelength = %f\n\tVias = %d\n\tMrd = %d\n",objective,total_cost,(int)total_vias,Mrd);
						
						
						// print the violated edges
						fprintf(outFile,"%d edge capacity violations\n",violations);
						for (i = 1; i <= m; i++)
							if ( (double)edge_count[i] > c[i] )
								fprintf(outFile,"\t%d : %d used, %d capacity\n",i,edge_count[i],c[i]);
						fprintf(outFile,"\n");
					
						// print the total cost of the solution
		//				fprintf(outFile,"Objective = %f\n\ttree cost = %f\n\tvias = %f\n",ALPHA*total_cost + BETA*total_vias,total_cost,total_vias);

						// print the actual solution
						print_all_nets(outFile,Solution,0,K);
						
						// close the output file
						fclose(outFile);
						
					} // fopen(outFile)
				}
				
			}// violations < best_violations
			
			tries++;
		
		}// violations > 0	
	} // multiple round
	
	
	free(x);


return 0;
}
