/*----------------------------------------------------------------*/
/* binheap.c                                                      */
/*    Implementation of a binary heap.  We define a structure that*/
/*  contains the key value, index, as well as array position.     */
/*  the position is contained so that we can locate in O(1) time. */
/*                                                                */
/*    We can build the heap by repeated inserts, or by feeding an */
/*  array of doubles to binheap_build.                            */
/*                                                                */
/*----------------------------------------------------------------*/

#ifndef BINHEAP_C_
#define BINHEAP_C_

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "binheap.h"
#include "fatal.h"

#define PARENT(i) i/2
#define LEFT(i) i*2 
#define RIGHT(i) i*2 + 1

/*-------------------------------------------------------------*/
/*	Purpose:  Initialize en empty heap                         */
/*  Returns:  Heap*                                            */
/*  Args   :                                                   */
/*		(int) MAX : Max num elements                           */
/*-------------------------------------------------------------*/
Heap *binheap_init(int MAX) {
	Heap *ret;
	
	ret = (Heap*) malloc(sizeof(Heap));
	ret->weight = (double*) malloc(sizeof(double)*(MAX+1));
	ret->index = (int*) malloc(sizeof(int)*(MAX+1));
	ret->position = (int*) malloc(sizeof(int)*(MAX+1));
	
	if ( (ret==NULL)||(ret->weight==NULL)||(ret->index==NULL)||(ret->position==NULL) )
		Error("binheap_init : memory allocation error\n");;
	
	return ret;
}


/*-------------------------------------------------------------*/
/*	Purpose: De-allocate a heap                                */
/*  Returns: nothing                                           */
/*  Args   :                                                   */
/*		(Heap*) h : pointer to the heap to be destroyed        */
/*-------------------------------------------------------------*/
void binheap_destroy(Heap *h) {
	free(h->weight);
	free(h->position);
	free(h->index);
	free(h);
}

/*-------------------------------------------------------------*/
/*	Purpose: (PRIVATE) exchange elements i and j in the heap   */                                                  
/*  Returns: nothing										   */
/*  Args   :												   */
/*		(Heap*) h : pointer to the heap						   */
/*		(int) i,j : elements to be swapped					   */
/*-------------------------------------------------------------*/
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

/*-------------------------------------------------------------*/
/*	Purpose: (PRIVATE) preserves min heap property of h		   */
/*  Returns: nothing										   */
/*  Args   :												   */
/*		(heap*) h ; pointer to the heap						   */
/*		(int)   i : position to heapify from				   */
/*-------------------------------------------------------------*/
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


/*-------------------------------------------------------------*/
/*	Purpose: Remove the smallest element from the heap.        */
/*  Returns: Nothing										   */
/*  Args   :												   */
/*		(Heap*) h : pointer to the heap						   */
/*-------------------------------------------------------------*/
void binheap_extract(Heap *h) {
	
	
	if (binheap_size(h) < 1) {
		printf("HEAP UNDERFLOW\n");
		exit(1);
	}

	binheap_exchange(h,1,binheap_size(h));

	binheap_size(h) = binheap_size(h) - 1;
	
	binheap_min_heapify(h,1);
}


/*-------------------------------------------------------------*/
/*	Purpose:  Build a min heap from an array of doubles        */
/*  Returns:  Heap* h (as argument)							   */
/*  Args   :												   */
/*		(Heap*) h   : pointer to the heap constructed		   */
/*		(double*) a : pointer to the array of key values       */
/*		(int) size  : size of array a						   */
/*	Other Notes	:											   */
/*		Must call binheap_init before calling build			   */
/*-------------------------------------------------------------*/
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

/*-------------------------------------------------------------*/
/*	Purpose: Insert an element into the heap				   */
/*  Returns: Nothing									       */
/*  Args   :												   */
/*		(Heap*) h	 : pointer to the heap					   */
/*		(double) key : key value of the new element			   */
/*-------------------------------------------------------------*/
int binheap_insert(Heap *h, double key) {
	
	binheap_size(h) = binheap_size(h) + 1;
	h->weight[binheap_size(h)] = DBL_MAX;
	h->index[binheap_size(h)] = binheap_size(h);
	h->position[binheap_size(h)] = binheap_size(h);
	
	binheap_decrease_key(h,binheap_size(h),key);
	
	return(binheap_size(h));
	
}

/*-------------------------------------------------------------*/
/*	Purpose: Display a heap to stdout						   */
/*  Returns: nothing										   */
/*  Args   :												   */
/*		(Heap*) h : pointer to the heap						   */
/*-------------------------------------------------------------*/
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
#endif
