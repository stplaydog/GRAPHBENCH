#include "graph.h"

struct elem{
	lint v_from;
	lint v_to;
	struct elem *next;
};

typedef struct elem *p_elem;
typedef struct elem t_elem;

void init_hash(lint num, p_elem *hash_map);
void insert_hash(lint from, lint to, p_elem *hash_map);
int search_hash(lint from, lint to, p_elem *hash_map);
