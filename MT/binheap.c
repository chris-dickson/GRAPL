#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define PARENT(i) i/2
#define LEFT(i) i*2 
#define RIGHT(i) i*2 + 1
#define binheap_indexofmin(h) h->index[1]
#define binheap_min(h) h->weight[1]
#define binheap_size(h) h->num_elements

typedef struct Heap_{
double	*weight;
int		*index;
int		*position;
int		num_elements;
}Heap;

/* IMPLEMENTATION */
Heap *binheap_init(int MAX) {
	Heap *ret;
	
	ret = (Heap*) malloc(sizeof(Heap));
	ret->weight = (double*) malloc(sizeof(double)*(MAX+1));
	ret->index = (int*) malloc(sizeof(int)*(MAX+1));
	ret->position = (int*) malloc(sizeof(int)*(MAX+1));
	
	return ret;
}

void binheap_destroy(Heap *h) {
	free(h->weight);
	free(h->position);
	free(h->index);
	free(h);
}

void binheap_exchange(Heap *h, int i, int j) {
	double	tmp_weight;
	int		tmp_index,tmp_pos;
	
	tmp_weight = h->weight[i];
	h->weight[i] = h->weight[j];
	h->weight[j] = tmp_weight;
	
	tmp_index = h->index[i];
	h->index[i] = h->index[j];
	h->index[j] = tmp_index;
	
	tmp_pos = h->position[ h->index[i] ];
	h->position[h->index[i]] = h->position[h->index[j]];
	h->position[h->index[j]] = tmp_pos;
}

void binheap_min_heapify(Heap *h, int i) {
	int l,r,smallest;
	
	l = LEFT(i);
	r = RIGHT(i);
	
	if (( l <= binheap_size(h) ) && ( h->weight[l] < h->weight[i] ))
		smallest = l;
	else
		smallest = i;
	
	if (( r <= binheap_size(h) ) && ( h->weight[r] < h->weight[smallest] ))
		smallest = r;
	
	if (smallest != i) {
		/* exchange */
		binheap_exchange(h,i,smallest);
		
		/* recurse */
		binheap_min_heapify(h,smallest);
	}
}

void binheap_extract(Heap *h) {
	
	
	if (binheap_size(h) < 1) {
		printf("HEAP UNDERFLOW\n");
		exit(1);
	}

	binheap_exchange(h,1,binheap_size(h));
	
/*	h->weight[1] = h->weight[ binheap_size(h) ]; */
	
	binheap_size(h) = binheap_size(h) - 1;
	
	binheap_min_heapify(h,1);
}

void binheap_build(Heap *h, double *A, int size) {
	int i;
	
	for (i = 1; i <= size; i++) {
		h->weight[i] = A[i];
		h->position[i] = i;
		h->index[i] = i;
	}
	
	binheap_size(h) = size;
	
	for (i = size/2; i >= 1; i--)
		binheap_min_heapify(h,i);
}

void binheap_decrease_key(Heap *h, int index, double key) {
	int i;
	
	i = h->position[index];
	
	
	h->weight[i] = key;
	while (( i > 1 ) && ( h->weight[PARENT(i)] > h->weight[i])) {
		binheap_exchange(h,i,PARENT(i));
		i = PARENT(i);
	}
}

int binheap_insert(Heap *h, double key) {
	
	binheap_size(h) = binheap_size(h) + 1;
	h->weight[binheap_size(h)] = DBL_MAX;
	h->index[binheap_size(h)] = binheap_size(h);
	h->position[binheap_size(h)] = binheap_size(h);
	
	binheap_decrease_key(h,binheap_size(h),key);
	
	return(binheap_size(h));
	
}

		
void binheap_print(Heap *h) {
	
	int i;
	printf("Displaying heap..\n");
	for (i = 1; i <= binheap_size(h); i++) 
		printf("%d : weight = %f , index = %d\n",i,h->weight[i],h->index[i]);
	
	printf("\nPositions:\n");
	for (i = 1; i <= binheap_size(h); i++)
		printf("Node %d in position %d\n",i,h->position[i]);
	printf("\n\n");
}
