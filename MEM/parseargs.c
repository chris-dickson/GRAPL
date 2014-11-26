/*--------------------------------------------------------------------------------------*/
/* parseargs.c																			*/
/*	Implementation of a simple parser to read the parameter file specified as an		*/
/* argument when the program is called.													*/
/*	To add/remove a parameter we have to change it in three places:						*/
/*		1.   This file : must add a symbol to the table, and add the if statement to	*/
/*			parse it.																	*/
/*		2.   parseargs.h : must add the parameter value, and parameter set variables	*/
/*			to the structure containing all the parameter settings						*/
/*		3.   main.c : must add a default value, and the if statement to parse set the   */
/*			actual value.  (Don't forget to create a varaible for it too)				*/
/*																						*/
/*																						*/
/*--------------------------------------------------------------------------------------*/


#ifndef PARSEARGS_C_
#define PARSEARGS_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "parseargs.h"


/*--------------------------------------------------------------*/
/*	Purpose: parse arguments from input file					*/
/*  Returns: (args*) ret										*/
/*  Args   :												    */
/*		(char*) arg_file : string holding name of argument file */
/*		(args*) ret		 : argument structure					*/
/*--------------------------------------------------------------*/

int parse_args(char *arg_file, args *ret) {
	FILE *fp;
	char line[LINE_LENGTH];
	
	/* symbol table*/
	char tau_sym[4]		="tau";
	char alpha_sym[6]	="alpha";
	char beta_sym[5]	="beta";
	char g_sym[2]		="g";
	char n_sym[2]		="n";
	char eps_sym[4]		="eps";
	char tauheur_sym[14]="tau_heuristic";
	char multiple_round_sym[16] ="multiple_round";
	char bail_on_lambda_sym[16] ="bail_on_lambda";
	char delta_lam_sym[10] ="delta_lam";
	char bench_style_sym[13] ="bench_style";
	char congestion_map_sym[15] ="congestion_map";
	char no_threads_sym[12]="no_threads";
	char warm_start_sym[12]="warm_start";

	if (( fp = fopen(arg_file,"r") ) == NULL) {
		printf("Error :  Cannot open file %s of parameter settings\nExiting.\n",arg_file);
		exit(1);
	}
	
	
	ret->set_N = 0;
	ret->set_TAU = 0;
	ret->set_ALPHA = 0;
	ret->set_BETA = 0;
	ret->set_TAU_HEURISTIC = 0;
	ret->set_g = 0;
	ret->set_EPS = 0;
	ret->set_MULTIPLE_ROUND = 0;
	ret->set_BAIL_ON_LAMBDA = 0;
	ret->set_DELTA_LAM = 0;
	ret->set_BENCH_STYLE = 0;
	ret->set_CONGESTION_MAP = 0;
	ret->set_NO_THREADS = 0;
	ret->set_WARM_START = 0;
	
	
	
	while (fgets(line,LINE_LENGTH,fp) != NULL) {
		char	c;
		char	var_name[LINE_LENGTH],var_value[LINE_LENGTH];
		int		i,j;
		int		found;
			
		for (i = 0; i < LINE_LENGTH; i++) {
			var_value[i] = '\0';
			var_name[i] = '\0';
		}
		
		
		i = 0;
		while ((c = line[i]) != '=') {
			var_name[i] = tolower(c);
			i++;
		}
		i++;			/* skip '=' symbol*/
		j = 0;
		while ((c = line[i]) != ';') {
			var_value[j] = tolower(c);
			j++;
			i++;
		}
		
		found = 0;
		
		/* Parse the current line*/
		if (strncmp(var_name,n_sym,LINE_LENGTH) == 0) {			/* N*/
			ret->N = atoi(var_value);
			ret->set_N = 1;
			found = 1;
		}
		if (strncmp(var_name,tau_sym,LINE_LENGTH) == 0) {		/* TAU*/
			ret->TAU = atof(var_value);
			ret->set_TAU = 1;
			found = 1;
		}		
		if(strncmp(var_name,alpha_sym,LINE_LENGTH) == 0) {		/* ALPHA*/
			ret->ALPHA = atof(var_value);
			ret->set_ALPHA = 1;
			found = 1;
		}
		if(strncmp(var_name,beta_sym,LINE_LENGTH) == 0) {		/* BETA*/
			ret->BETA = atof(var_value);
			ret->set_BETA = 1;
			found = 1;
		}
		if(strncmp(var_name,g_sym,LINE_LENGTH) == 0) {			/* g*/
			ret->g = atof(var_value);
			ret->set_g = 1;
			found = 1;
		}
		if(strncmp(var_name,tauheur_sym,LINE_LENGTH) == 0) {	/* TAU_HEURISTIC*/
			ret->TAU_HEURISTIC = atoi(var_value);
			ret->set_TAU_HEURISTIC = 1;
			found = 1;
		}
		if (strncmp(var_name,eps_sym,LINE_LENGTH) == 0) {		/* EPS*/
			ret->EPS = atof(var_value);
			ret->set_EPS = 1;
			found = 1;
		}
		if (strncmp(var_name,multiple_round_sym,LINE_LENGTH) == 0) {	/*MULTIPLE_ROUND*/
			ret->MULTIPLE_ROUND = atoi(var_value);
			ret->set_MULTIPLE_ROUND = 1;
			found = 1;
		}
		if (strncmp(var_name,bail_on_lambda_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			ret->BAIL_ON_LAMBDA = atof(var_value);
			ret->set_BAIL_ON_LAMBDA = 1;
			found = 1;
		}
		if (strncmp(var_name,delta_lam_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			ret->DELTA_LAM = atof(var_value);
			ret->set_DELTA_LAM = 1;
			found = 1;
		}
		if (strncmp(var_name,bench_style_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			ret->BENCH_STYLE = atoi(var_value);
			ret->set_BENCH_STYLE = 1;
			found = 1;
		}

		if (strncmp(var_name,congestion_map_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			ret->CONGESTION_MAP = atoi(var_value);
			ret->set_CONGESTION_MAP = 1;
			found = 1;
		}

		if (strncmp(var_name,no_threads_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			ret->NO_THREADS = atoi(var_value);
			ret->set_NO_THREADS = 1;
			found = 1;
		}

		if (strncmp(var_name,warm_start_sym,LINE_LENGTH) == 0) {	/*BAIL_ON_LAMDA*/
			strcpy(ret->WARM_START,var_value);
			ret->set_WARM_START = 1;
			found = 1;
		}


		if (!found) {
			printf("Error parsing argument file.\nUnknown symbol:\n\t %s \n",var_name);
			exit(1);
		}
		
	}
	fclose(fp);
	return 0;
}

#endif
