#include <stdio.h>
#include <stdlib.h>

#define N 10
#define MAX_INT 3000

static int qsort_match(const void *a1, const void *a2) {
	if (*(int*)a1 < *(int*)a2)
		return -1;
		else if (*(int*)a1 > *(int*)a2)
			return 1;
	return 0;
}

int main() {
	int i;
	int *SteinerTree;
	
	SteinerTree = (int*) malloc(sizeof(int)*N);
	
	for (i = 0; i < N; i++) 
		SteinerTree[i] = (int) ( rand() % MAX_INT ) + 1;
		
	for (i = 0; i < N; i++)
		printf("%d ",SteinerTree[i]);
	printf("\n");

	qsort(SteinerTree , N, sizeof(int), &qsort_match);
	
		for (i = 0; i < N; i++)
		printf("%d ",SteinerTree[i]);
	printf("\n");

	
	
	return 0;
}
