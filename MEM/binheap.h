/*----------------------------------------------------------------*/
/* binheap.h                                                      */
/*     Header file for a binary heap.                             */
/*----------------------------------------------------------------*/

#ifndef BINHEAP_H_
#define BINHEAP_H_

typedef struct Heap_{
	double	*weight;
	int		*index;
	int		*position;
	int		num_elements;
}Heap;

Heap *binheap_init(int MAX);
void binheap_build(Heap *h, double *A, int size);
void binheap_destroy(Heap *h);
void binheap_extract(Heap *h);
void binheap_decrease_key(Heap *h, int index, double key);
int binheap_insert(Heap *h, double key);
void binheap_print(Heap *h);
#define binheap_indexofmin(h) h->index[1]
#define binheap_min(h) h->weight[1]
#define binheap_size(h) h->num_elements
#define binheap_isempty(h) (h->num_elements == 0)

#endif
