#include "hashmap.h"
#include "utils.h"


void init_hash(lint num, p_elem *hash_map){
	lint i;
	for (i=0;i<num;i++){
		hash_map[i] = (p_elem) malloc(sizeof(t_elem));
		hash_map[i]->v_from =0;
		hash_map[i]->v_to =0;
		hash_map[i]->next =NULL;
	}
}


void insert_hash(lint from, lint to, p_elem *hash_map){
	int hash_val = hash(from+to);
	p_elem e = (p_elem)malloc(sizeof(t_elem));
	e->v_from =from;
	e->v_to = to;
	e->next = NULL;
	if(hash_map[hash_val]->next==NULL){
		hash_map[hash_val]->next = e;	
	}	
	else{
		p_elem tmp  = hash_map[hash_val]->next;
		while(tmp->next !=NULL)
			tmp = tmp->next;
		tmp->next = e;
	}
}
int search_hash(lint from, lint to, p_elem *hash_map){
	int hash_val = hash(from+to);
	if(hash_map[hash_val]->next==NULL){
		return FALSE;
	}	
	else{
		p_elem tmp = hash_map[hash_val]->next;
		while(TRUE){
			if(tmp->v_from==from && tmp->v_to == to){
				//DPRINTF(1, "DUPLICATE %lld %lld\n", from, to);
				return TRUE;
			}
			else if(tmp->next != NULL)
				tmp = tmp->next;
			else
				break;
		}
	}
	return FALSE;
}

void free_hash(lint num, p_elem *hash_map){
	lint i;
	for(i=0;i<num; i++)
		free(hash_map[i]);
	free(hash_map);
}
