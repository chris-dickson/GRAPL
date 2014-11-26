/*--------------------------------------------------------------------------------------*/
/* routing.h																			*/
/*		This is the header for some standard routines for our routing problem.  It		*/
/* contains debug routines, as well as routines to parse the benchmarks.				*/
/*																						*/
/*		This file also defines the NET data structure.  This is the data type for our	*/
/*	main data for routing																*/
/*--------------------------------------------------------------------------------------*/


#ifndef ROUTING_H_
#define ROUTING_H_

#include "helper.h"
#include "set.h"
#include "graph.h"

#define IBM_BENCH 0
#define YAL_BENCH 1

#define LINELENGTH 1024
#define ID_LENGTH 10

#define MAX_NO_EDGES 2048


typedef struct Net_ {
	char		ID[ID_LENGTH];
	int			net_num;
	int			no_terminals;
	CoordData	*terminals;
	int			*SteinerTree;
	int			edge_count;
	int			vias;
	double		tree_cost;
	int			routed;		/* boolean, 0 if not routed*/
	int			fixed;
}Net;

int parseBenchmark(char *fName, int *no_nets, Net **T, int N, int *gw, int *gh, int *hCap, int *vCap, int *num_routed);
void parse_yal_benchmark(char *fName, int *no_nets, Net **T, int N, int *gw, int *gh, int *hCap, int *vCap, double *hLength, double *vLength);


int U(Net *T, int edge_i, int tree_j);
void find_unique_terminals(Net *T, int no_nets, CoordData **unique_terminals, int *no_unique_terminals);
/* debug*/
void print_all_nets(FILE *fp,Net *T, int N, int K);
void print_one_net(FILE *fp, Net *T, int k);
void print_range_net(FILE *fp, Net *T, int lb, int ub);
void dump_Array(double *a, const char *prefix, int extension, int size);
void query_mode(Net *T,int K);
void edge_Histogram(Graph *G, double *f, int base, int size, const char *prefix, int extension);


void fix_net(Net *T, int i);

void get_next_tree(FILE *fp, Graph *g, double hl, double vl, int **SteinerTree, int *edge_count_OUT, double *tree_cost_OUT, char *ID_OUT);


#endif
