/*----------------------------------------------------------------------------------*/
/* frac_round_MRD.c																	*/
/*		This is the main function to the rounding program.  This program will round */
/*	the [0,1]-valued solution to a {0,1} (integer) solution.						*/
/*																					*/
/*		This program accepts three command line arguments.   The first argument is	*/
/*	the file containing the trees generated.  The second argument is the stat file	*/
/*	produced.   The third argument is the solution file that is to be written.		*/
/*	This program will parse the stat file to obtain the step length and generate	*/
/*	the corresponding indicator (x) vector.  It will then read in the trees and		*/
/*	begin rounding.																	*/
/*																					*/
/*		This version of the rounding code will minimize the maximum routing density */
/*	(MRD, otherwise known as the maximum edge congestion).  The program always keeps*/
/*	the best solution written to the output file.  We judge the 'best' solution by	*/
/*	the following scheme:															*/
/*		while (not done rounding) do												*/
/*			let x be the solution obtained from randomized rounding					*/
/*			if mrd(x) < mrd(best) and best is infeaseable then						*/
/*				best = x															*/
/*			else if best is feaseable and x is feaseable then						*/
/*				if objective(x) < objective(best) then								*/
/*					best = x														*/
/*		end																			*/
/*																					*/
/*		In plain words, we always try to minimize the MRD until we get a feaseable	*/
/*	solution.   Once we get a feaseable solution, we try to minimize the objective. */
/*----------------------------------------------------------------------------------*/


#define MAX_TRIES 10000

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "list.h"
#include "set.h"
#include "graph.h"
#include "graphalg.h"
#include "rounding.h"


#define MAX_WORD_SIZE 100
#define LINE_SIZE 1000	
#define MULTIPLE_ROUND 1
#define ID_LENGTH 256

typedef struct Net_ {
	char		ID[10];
	int			net_num;
	int			no_terminals;
	CoordData	*terminals;
	int			*SteinerTree;
	int			edge_count;
	int			vias;
	double		tree_cost;
	int			routed;		/* boolean, 0 if not routed */
}Net;

#include "frac_round_aux.c"

int main(int argc, char *argv[]) {
	FILE	*inSolutionFile, *inStatFile, *outFile;
	char	a_word[MAX_WORD_SIZE],a_line[LINE_SIZE];
	int		i,done,iters,l,k;
	
	int		Mrd;
	int		best_Mrd;
	
	int		N,K;
	double	*x;
	int		*x_star;
	Net		*T,*Solution;
	
	Graph	G,H;
	int		width,height;
	int		hc,vc;
	int		*c;
	int		m;
	
	double	ALPHA,BETA;
	char	float_value[32];
	double	total_cost,total_vias;
	double	objective;
	double	best_objective;
	
	time_t		iteration_start,iteration_stop;
	double		iteration_ET;
	int			hours,mins,sec;
	
	int			time_limit;
	
	
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

	printf("How many iterations?\n");
	scanf("%d",&iters);
	
	printf("Enter time limit (0 for unlimited): ");
	scanf("%d",&time_limit);
	if (time_limit < 0)
		fprintf(stderr,"time limit must be non-negative\nprogram exiting...\n"),exit(1);

	/* figure out width */
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Width")) {
			done = 1;
		}
	}
	fscanf(inStatFile,"%d",&width);
	
	/* figure out height */
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Height")) {
			done = 1;
		}
	}
	fscanf(inStatFile,"%d",&height);
	
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Vertical")) {
			
			fscanf(inStatFile,"%s",a_word);
			if (!strcmp(a_word,"Capacity"))
				done = 1;
		}
	}
	fscanf(inStatFile,"%d",&vc);

	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Horizontal")) {
			
			fscanf(inStatFile,"%s",a_word);
			if (!strcmp(a_word,"Capacity"))
				done = 1;
		}
	}
	fscanf(inStatFile,"%d",&hc);
	
	printf("Alpha: ");
	scanf("%s",float_value);
	ALPHA = atof(float_value);
	
	printf("Beta: ");
	scanf("%s",float_value);
	BETA = atof(float_value);
	
	srand( time(NULL) );
	

	/* figure out K */
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Num")) {
			done = 1;
		}
	}
	fscanf(inStatFile,"%*s %d",&K);
	
	/* figure out N */
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
	

	/* find first line */
	done = 0;
	while (!done) {
		fscanf(inStatFile,"%s",a_word);
		if (!strcmp(a_word,"Iteration")) {
			done = 1;
		}
	}
	fgets(a_line,LINE_SIZE,inStatFile); /* rest of first line */
	fgets(a_line,LINE_SIZE,inStatFile); /* second line, TAU = 1.0 */
	
	for (i = 1; i <= K; i++)
		x[i] = 1.0;

	done = 0;
	while (!done) {
		char *ret;
		
		ret = fgets(a_line,LINE_SIZE,inStatFile);
		if (ret == NULL)							/* END OF FILE */
			done = 1;
		else {
			int iteration;
			char TAU[MAX_WORD_SIZE];
			double TAU_VAL;
			int		k;
			
			sscanf(a_line,"%d %*s %*s %s %*s %*s %*s %*s %*s %*s\n",&iteration,TAU);
			TAU_VAL = atof(TAU);
			
			for (l = 0; l < iteration; l++)
				for (k = 1; k <= K; k++)
					x[l*K + k] = x[l*K + k]*(1-TAU_VAL);
			
			for (k = 1; k <= K; k++)
				x[iteration*K + k] = TAU_VAL;
			
			if (iteration == iters)
				done = 1;
		}
	}
	
	/* We now have our x array.   Now we must parse the solution file before we attempt to round */
	if ((T = (Net*) malloc(sizeof(Net)*((iters+1)*K + 1)))==NULL) {
		printf("benchmark, temp : memory allocation error\n"); fflush(stdout);
		exit(1);
	}

		
	
	/* parse the trees */
	for (l = 0; l <= iters; l++) {
		int iter;
		
		fscanf(inSolutionFile,"%*s %d",&iter);		/* Iteration %d			*/
		fscanf(inSolutionFile,"%*s");				/* -------------------	*/
		fscanf(inSolutionFile,"%*s");				/* \n					*/
		
		/* for each tree	*/
		for (k = 1; k <= K; k++) {
			char no_term_string[MAX_WORD_SIZE], coord[MAX_WORD_SIZE], ID[ID_LENGTH];
			
			fscanf(inSolutionFile,"%*s %*d");							/* Net %d				*/
		    fscanf(inSolutionFile,"%*s %s",ID);							/* ID %s				*/
			strcpy(T[l*K + k].ID,ID);
			fscanf(inSolutionFile,"%s",no_term_string);					/* no_terminals=%d		*/
			i = 0;														/* skip to the = sign	*/
			while (no_term_string[i] != '=')
				i++;
			i++;														/* point to the number	*/
			T[l*K+k].no_terminals = atoi( &(no_term_string[i]) );		/* convert				*/
	
			T[l*K+k].terminals = (CoordData*) malloc(sizeof(CoordData)*T[l*K+k].no_terminals);
				   

			fscanf(inSolutionFile,"%*s %*s");							/* terminal set			*/
			for (i = 0; i < T[l*K+k].no_terminals; i++) {
				
				fscanf(inSolutionFile,"%*d %*s %s",coord);				/* %d : (%d,%d,%d)		*/
				sscanf(coord,"(%hd,%hd,%hd)",&(T[l*K+k].terminals[i].x),&(T[l*K+k].terminals[i].y),&(T[l*K+k].terminals[i].z));				
			}
			
			/* get edge count	*/
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].edge_count = atoi(&a_line[i]);
			
			/* get vias			*/
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].vias = atoi(&a_line[i]);
			
			/* get tree cost	*/
			fscanf(inSolutionFile,"%s",a_line);
			i = 0;
			while (a_line[i] != '=')
				i++;
			i++;
			T[l*K+k].tree_cost = atof(&a_line[i]);
			
			fgets(a_line,LINE_SIZE,inSolutionFile);		/*	junk		*/
			fgets(a_line,LINE_SIZE,inSolutionFile);		/*	steiner tree*/
			
			T[l*K+k].SteinerTree = (int*) malloc(sizeof(int)*T[l*K+k].edge_count);
			

			for (i = 0; i < T[l*K+k].edge_count; i++) 
				fscanf(inSolutionFile,"%d",&(T[l*K+k].SteinerTree[i]));
			
			qsort(T[l*K+k].SteinerTree , T[l*K+k].edge_count , sizeof(int), &qsort_match);
				

			
			fgets(a_line,LINE_SIZE,inSolutionFile);			
			fgets(a_line,LINE_SIZE,inSolutionFile);			
			
			
		}
	}
	
			
	/* We now have our fractional solution stored in T as it would be in route.   Now, we must make the graphs	*/
	fclose(inSolutionFile);
	fclose(inStatFile);
	
	gen_gridgraph(&G,width,height);
	grid2VL(&G,&H,width,height);
	if ((c = (int*) malloc(sizeof(int)*(((G.ecount)/2)+1+1)))==NULL) {
		printf("c : memory allocation error\n"); fflush(stdout);
		exit(1);
	}
	set_capacity_lookup(&c,&H,hc,vc);

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
		
		best_objective = DBL_MAX;
		best_violations = m;
		violations = m;
		best_Mrd = m;
		
		tries = 1;
		iteration_start = time(NULL);
		while (tries < MAX_TRIES) {
			
			/* perform a randomized rounding	*/
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
			
			/* compute congestion	*/
			for (k = 1; k <= K; k++) 
				for (i = 1; i <= m; i++) 
					edge_count[i] += U_sorted(T,i,x_star_inverse[k]);
			
			/* check for edge violations	*/
			Mrd = 0;
			for (i = 1; i <= m; i++) {
				if ((double)edge_count[i] > c[i])
					violations++;
				if (edge_count[i] > Mrd)
					Mrd = edge_count[i];
			}
			
			/* compute objective value		*/
			total_cost = 0.0;
			total_vias = 0.0;
			for (k = 1; k <= K; k++) {
				total_cost += (double) T[x_star_inverse[k]].tree_cost;
				total_vias += (double) T[x_star_inverse[k]].vias;
			}
			objective = ALPHA*total_cost + BETA*total_vias;			
			
			iteration_stop = time(NULL);
			
			/* is this round any better?	*/
			if ((Mrd < best_Mrd)||((Mrd == best_Mrd)&&(objective < best_objective))) {

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
					printf("\tObjective = %E (After %d tries and %dh %dm %ds, violations = %d)\n",objective,tries,hours,mins,sec,violations);
				

					if ((outFile = fopen(argv[3],"w"))!=NULL) {
						
						/* overwrite old solution	*/
						for (k = 1; k <= K; k++)
							memcpy( &Solution[k] , &T[x_star_inverse[k]] , sizeof(T[x_star_inverse[k]]) );
						
						/* print header				*/
						if (violations > 0) fprintf(outFile,"INFEASEABLE INTEGER SOLUTION\n");
						else fprintf(outFile,"INTEGER SOLUTION\n");
						fprintf(outFile,"----------------------------\n");
						
						fprintf(outFile,"OBJECTIVE VALUE = %E\n\tWirelength = %f\n\tVias = %d\n\tMrd = %d\n",objective,total_cost,(int)total_vias,Mrd);
						fprintf(outFile,"%d hours, %d mins, %d seconds\n",hours,mins,sec);
						
						
						/* print the violated edges	*/
						fprintf(outFile,"%d edge capacity violations out of %d total edges (%lf percent)\n",violations,width*(height-1) + height*(width-1), ( (double)violations / (width*(height-1) + height*(width-1)))*100);
						for (i = 1; i <= m; i++)
							if ( (double)edge_count[i] > c[i] )
								fprintf(outFile,"\t%d : %d used, %d capacity\n",i,edge_count[i],c[i]);
						fprintf(outFile,"\n");
					
						/* print the actual solution	*/
						print_all_nets(outFile,Solution,0,K);
						
						/* close the output file		*/
						fclose(outFile);
						
					} /* fopen(outFile)	*/
				}
				
			}/* violations < best_violations	*/
	
			tries++;
			iteration_ET = difftime(iteration_stop,iteration_start);
			hours = ((int)iteration_ET) / 3600;
			mins  = ((int)iteration_ET - (3600*hours)) / 60;
			sec   = (int)iteration_ET - (3600*hours) - 60*mins;
			
			/* timed out, stop rounding			*/
			if ((mins >= time_limit)&&(time_limit>0)) {
				printf("TIME LIMIT REACHED\n");
				tries = MAX_TRIES;
			}
		
		}/* violations > 0		*/
	} /* multiple round		*/

	free(x);

return 0;
}
