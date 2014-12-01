/*--------------------------------------------------------------------------------------*/
/* prune_yal.c																			*/
/*		This program prunes yal style benchmarks.	It removes repeated terminals, as	*/
/*	well as allowing us to provide an upper bound on the number of terminals for each	*/
/*	net.																				*/
/*--------------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "set.h"
#include "list.h"

typedef struct CoordData_ {
	int x,y;
}CoordData;

typedef struct Net_ {
//  char *net_name;
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
	int		i,j,k;
	Net		*temp;
	char	tmp_str[64];
	
	//benchmark data
	int		gw,gh;
	double	hl,vl;
	int		vCap,hCap;
	int		no_nets;
	int		no_modules;
	
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
		
	
	fscanf(fp,"%*s %*s %*s %d",&no_nets);			
	fscanf(fp,"%*s %*s %*s %d",&no_modules);					
	fscanf(fp,"%*s %*s %*s %d",&gw);					
	fscanf(fp,"%*s %*s %*s %d",&gh);					
	fscanf(fp,"%*s %*s %*s %d",&hCap);			
	vCap = hCap;
	
	fscanf(fp,"%*s %*s %*s %s",tmp_str);
	hl = atof(tmp_str);
	
	fscanf(fp,"%*s %*s %*s %s",tmp_str);
	vl = atof(tmp_str);
	
	temp = (Net*) malloc(sizeof(Net)*no_nets);
	
	routed = 0;
	for (k = 0; k < no_nets; k++) {
		int i;
		Set	pointSet;
		int			j;
		int			tmp_net_num;
		ListElmt	*e;

		
		fscanf(fp,"%d, %d:,",&tmp_net_num,&(temp[k].no_terminals));
		
		if (temp[k].no_terminals == 0) {
			fscanf(fp,"%*s");
			continue;
		}
		
		set_init( &pointSet, match2d, free );

		
		for (i = 0; i < temp[k].no_terminals; i++) {
			char a_string[16];
			point2d *tmp;
			int		x,y;

			
			fscanf(fp,"%s",a_string);
			sscanf(a_string,"%d,%d,",&x,&y);
			
			if ((tmp = (point2d*) malloc(sizeof(point2d)))==NULL){
				printf("benchmarks.c : tmp mem allocation error\n");
				fflush(stdout);
				exit(1);
			}
			tmp->x = x-1;	tmp->y = y-1;
			
			set_insert(&pointSet,tmp);

			if (i == temp[k].no_terminals - 1) 
				fgets(a_string,16,fp);			// pick up the trashŒ
			
		}
		
		if ((temp[k].terminals = (CoordData*) malloc(sizeof(CoordData) * set_size(&pointSet))) == NULL) {
			printf("benchmark, temp[i].terminals : memory allocation error\n"); fflush(stdout);
			fflush(stdout);
			exit(1);
		}
		
		temp[k].no_terminals = set_size(&pointSet);
		temp[k].net_num = tmp_net_num;
		
		if ((temp[k].no_terminals > 1)&&(temp[k].no_terminals <= cut_off))
			routed++;

		j = 0;
		for (e = list_head(&pointSet); e != NULL; e = list_next(e)) {
			point2d		*data;
			
			data = (point2d*) list_data(e);
			temp[k].terminals[j].x = data->x; 
			temp[k].terminals[j].y = data->y;
			j++;
		}		
	}
	fclose(fp);
					
/*------------------------------------------------------------------------*/
/* Reconstruct benchmark, but with pruned entries						  */
/*------------------------------------------------------------------------*/
	fprintf(fpOut,"Number of Nets: %d\n",routed);
	fprintf(fpOut,"Number of Modules: %d\n",no_modules);
	fprintf(fpOut,"Number of Rows: %d\n",gw);
	fprintf(fpOut,"Number of Columns: %d\n",gh);
	fprintf(fpOut,"The Maximum Capacity: %d\n",hCap);
	fprintf(fpOut,"Horizontal Edge Length: %f\n",hl);
	fprintf(fpOut,"Vertical Edge Length: %f\n\n\n",vl);	
	
	for (i = 0; i < no_nets; i++) {
		if ((temp[i].no_terminals > 1)&&(temp[i].no_terminals <= cut_off)) {
	
			fprintf(fpOut,"%d, %d:, ",temp[i].net_num,temp[i].no_terminals);
			for (j = 0; j < temp[i].no_terminals; j++) 
				fprintf(fpOut,"%d,%d, ",temp[i].terminals[j].x,temp[i].terminals[j].y);
			fprintf(fpOut,";\n");
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
