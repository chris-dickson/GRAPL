#ifndef HASH_H_
#define HASH_H_

#include <stdio.h>
#include <stdlib.h>

typedef struct hash_node_{
	int x,y,z;
	struct hash_node_ *next;
} hash_node;

typedef hash_node hash_list;

typedef struct hash_table_{
	hash_list	**D;
	int			no_lists;
	
	int			(*h)(int x);
	
}hash_table;

hash_list*
hash_list_init() {
	return NULL;
}

hash_list*
hash_list_insert(hash_list *L, int x, int y, int z) {
	hash_list *n;
	
	n = (hash_list*) malloc(sizeof(hash_list));
	n->x = x;
	n->y = y;
	n->z = z;
	n->next = L;
	L = n;
	return L;
}

void
hash_table_init(hash_table *H, int N) {
	int i;
	
	H->D = (hash_list**) malloc(sizeof(hash_list*)*N);
	for (i = 0; i < N; i++)
		H->D[i] = hash_list_init();
}


	




#endif
