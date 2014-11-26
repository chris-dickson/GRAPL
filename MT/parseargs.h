/*--------------------------------------------------------------------------------------*/
/* parseargs.h																			*/
/*	Header file for a simple parser to read the parameter file specified as an			*/
/* argument when the program is called.													*/
/*																						*/
/* Chris Dickson																		*/
/* McMaster University																	*/
/* June, 2006																			*/
/*--------------------------------------------------------------------------------------*/



#ifndef PARSEARGS_H_
#define PARSEARGS_H_

#define LINE_LENGTH 1024

#define BB_AREA 1
#define BB_SUM 2

#define INCREASING 1
#define DECREASING 0


typedef struct args_{
	int		N;
	double	ALPHA,BETA,g,TAU;
	int		TAU_HEURISTIC;
	double	EPS;
	int		MULTIPLE_ROUND;			/* boolean, used if we keep rounding*/
	double	BAIL_ON_LAMBDA;			/* quits if lambda lower than this value*/
	double	DELTA_LAM;				/* quits if change in lambda is less than this value*/
	int		BENCH_STYLE;
	int		CONGESTION_MAP;
	int		NO_THREADS;
	char	WARM_START[LINE_LENGTH];
	char	STAT_OUTPUT[LINE_LENGTH];
	char	TREE_OUTPUT[LINE_LENGTH];
	int		R;
	int		FIXING_HEURISTIC;
	double	FIXING_PERCENTAGE;
	int		FIXING_ORDER;
	
	int set_N,set_ALPHA,set_BETA,set_g,set_TAU,set_TAU_HEURISTIC,set_EPS,set_MULTIPLE_ROUND,
		set_BAIL_ON_LAMBDA, set_DELTA_LAM, set_BENCH_STYLE, set_CONGESTION_MAP, set_NO_THREADS, set_WARM_START, set_TREE_OUTPUT, 
		set_STAT_OUTPUT, set_FIXING_HEURISTIC, set_FIXING_PERCENTAGE, set_R, set_FIXING_ORDER;
	
}args;

int parse_args(char *arg_file, args *ret);

#endif
