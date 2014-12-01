/*--------------------------------------------------------------------------------------*/
/* prune_ibm.c																			*/
/*		This program prunes ibm style benchmarks.	It removes repeated terminals, as	*/
/*	well as allowing us to provide an upper bound on the number of terminals for each	*/
/*	net.																				*/
/*--------------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "set.h"
#include "list.h"

typedef struct CoordData_ {
	int x,y;
}CoordData;

typedef struct Net_ {
  char net_name[64];
	int net_num;
	int no_terminals;
	CoordData* terminals;
}Net;

typedef struct point2d_ {
	int x,y;
} point2d;

int match2d(const void *a1, const void *a2) {
	return ((((point2d*)a1)->x == ((point2d*)a2)->x)&&(((point2d*)a1)->y == ((point2d*)a2)->y));
}

/*--------------------------------------------------------------*/
/*	Purpose: Remove all repeated terminals from each net.  Also	*/
/*			we can include an upper bound on the number of		*/
/*			terminals in each net								*/
/*  Returns: output file, specified as second argument			*/
/*  Args   :												    */
/*		(char*) inFileName,outFileName  : input/output files	*/
/*		(int) Tolerance : upper bound on # of terminals / net	*/
/*--------------------------------------------------------------*/
int prune_benchmark(char *inFileName, char *outFileName, int Tolerance) {
	FILE	*fp,*fpOut;
	int		i,j;
	Net		*temp;
	
	//benchmark data
	int		gw,gh;
	int		vCap,hCap;
	int		no_nets;
	char	net_name[64];
	
	int		routed;
	int		cut_off;
	
	if (Tolerance == 0)
		cut_off = INT_MAX;
	else
		cut_off = Tolerance;
	
	


	// open the benchmark file
	if ( (fp = fopen(inFileName,"r")) == NULL) {
		printf("Can't find the file "); printf(inFileName); printf("for reading\n");
		return 1;
	}
	
	if ( (fpOut = fopen(outFileName,"w")) == NULL) {
		printf("Can't open the file "); printf(outFileName); printf(" for writing\n");
		return 1;
	}
		
	
	// grid dimensions
	fscanf(fp,"%*s %d %d",&gw,&gh);
	
	// vertical capacity
	fscanf(fp,"%*s %*s %d", &vCap);
	
	// horizontal capacity
	fscanf(fp,"%*s %*s %d", &hCap);
	
	// num nets
	fscanf(fp,"%*s %*s %d",&no_nets);
	
	temp = (Net*) malloc(sizeof(Net)*no_nets);
			
	routed = 0;
	
	//begin parsing nets
	printf("Parsing nets...");
	for (i = 0; i < no_nets; i++) {
		Set			pointSet;
		int			numPoints,j;
		int			tmp_net_num;
		ListElmt	*e;
		
		// initalize a point set
		set_init( &pointSet, match2d, free );
		
	//	printf("Parsing net %d\n",i);
		fscanf(fp,"%s %d %d",net_name,&tmp_net_num,&numPoints);
		
		
		// throw all points into a set, duplicates will get removed
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
			
			set_insert(&pointSet,tmp);
		}
		
		
		// convert point set to an array of CoordData
		if ((temp[i].terminals = (CoordData*) malloc(sizeof(CoordData) * set_size(&pointSet))) == NULL) {
			printf("benchmark, temp[i].terminals : memory allocation error\n"); fflush(stdout);
			fflush(stdout);
			exit(1);
		}
		
		strcpy(temp[i].net_name,net_name);
		temp[i].no_terminals = set_size(&pointSet);
		temp[i].net_num = tmp_net_num;
		
		if ((temp[i].no_terminals > 1)&&(temp[i].no_terminals <= cut_off))
			routed++;

		j = 0;
		for (e = list_head(&pointSet); e != NULL; e = list_next(e)) {
			point2d		*data;
			
			data = (point2d*) list_data(e);
			temp[i].terminals[j].x = data->x; 
			temp[i].terminals[j].y = data->y;
			j++;
		}		
	}
	printf("done\n");
	fclose(fp);
	
	
/* Reconstruct benchmark, but with pruned entries*/
	fprintf(fpOut,"grid %d %d\n",gw,gh);
	fprintf(fpOut,"vertical capacity %d\n",vCap);
	fprintf(fpOut,"horizontal capacity %d\n",hCap);
	fprintf(fpOut,"num net %d\n",routed);
	
	
	for (i = 0; i < no_nets; i++) {
		if ((temp[i].no_terminals > 1)&&(temp[i].no_terminals <= cut_off)) {
	
			
			fprintf(fpOut,"%s %d %d \n",temp[i].net_name,temp[i].net_num,temp[i].no_terminals);
			for (j = 0; j < temp[i].no_terminals; j++) 
				fprintf(fpOut,"  %d %d\n",temp[i].terminals[j].x,temp[i].terminals[j].y);
		}
	}
	fclose(fpOut);
		
	return 0;
}	


/*--------------------------------------------------------------*/
/*	Purpose: call pruning function based on command line args	*/
/*  Returns: 0 always											*/
/*  Usage   :	see below										*/
/*--------------------------------------------------------------*/
int main(int argc, char *argv[]) {
	
	if (argc == 3) 
		prune_benchmark(argv[1],argv[2],0);
	else if (argc == 4)
		prune_benchmark(argv[1],argv[2],atoi(argv[3]));
	else {
		printf("Usage :  prune_benchmark file1 file2 [tol]\n");
		printf("Where :\n");
		printf("\tfile1\t\tThe benchmark file to be pruned\n");
		printf("\tfile2\t\tThe pruned benchmark file to be written to\n");
		printf("\ttol\t\tOptional tolerance, only nets less than or equal to tol\n");
		printf("\t\t\twill be output to file2.   If left blank, all nets are output\n");
		return 1;
	}
	
	return 0;
}
