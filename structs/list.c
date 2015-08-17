#include "list.h"

void
creat_list(p_list list){
	list->int_val = -1;
	list->int_val1 = -1;
	list->double_val = -1;
	list->next = NULL;	
}

void
add_list_one(p_list list, lint val_int){
	p_list elem = (p_list)malloc(sizeof(t_list));
	elem->int_val = val_int;
	elem->int_val1 = -1;
	elem->double_val = -1;
	elem->next = NULL;	
	
	if(list->next==NULL){
		list->next = elem;	
	}else{
		p_list tmp = list->next;
		elem->next = tmp;
		list->next = elem;
	}
}

void
add_list_two(p_list list, lint val_int, double val_double){
	p_list elem = (p_list)malloc(sizeof(t_list));
	elem->int_val = val_int;
	elem->int_val1 = -1;
	elem->double_val = val_double;
	elem->next = NULL;	
	
	if(list->next==NULL){
		list->next = elem;	
	}else{
		p_list tmp = list->next;
		elem->next = tmp;
		list->next = elem;
	}
}

void
add_list_three(p_list list, lint val_int, lint val1_int, double val_double){
	p_list elem = (p_list)malloc(sizeof(t_list));
	elem->int_val = val_int;
	elem->int_val1 = val1_int;
	elem->double_val = val_double;
	elem->next = NULL;	
	
	if(list->next==NULL){
		list->next = elem;	
	}else{
		p_list tmp = list->next;
		elem->next = tmp;
		list->next = elem;
	}
}

p_list
pop_elem_list(p_list list){
	p_list tmp = list->next;
	if(tmp==NULL)
		return NULL;
	list->next = tmp->next;
	return tmp;
} 

p_list
pop_min(p_list list){
	p_list tmp = list->next;
	if(tmp==NULL)
		return NULL;
	double min = INFINITY;	
	lint min_idx = -1;
	while(tmp != NULL){
		double val = tmp->double_val;
		lint idx = tmp->int_val;
		if(val<min){
			min = val;
			min_idx = idx;	
		}
		tmp = tmp->next;
	}
	tmp = list->next;
	p_list cur = list;
	while(tmp != NULL){
		double val = tmp->double_val;
		lint idx = tmp->int_val;
		if(idx == min_idx){
			cur->next = tmp->next;	
			return tmp;
		}
		tmp = tmp->next;
		cur = cur->next;
	}	
	return NULL;
}

p_list
pop_find(p_list list, lint query){
	p_list tmp = list->next;
	p_list pre = list;
	if(tmp==NULL)
		return NULL;
	while(tmp != NULL){
		lint idx = tmp->int_val;
		if(idx==query){
			pre->next = tmp->next;
			return tmp;
		}
		tmp = tmp->next;
		pre = pre->next;
	}
	return NULL;
}

sint 
find(p_list list, lint query){
	p_list tmp = list->next;
	while(tmp != NULL){
		lint idx = tmp->int_val;
		if(idx == query)
			return TRUE;
		tmp = tmp->next;
	}
	return FALSE;
}

void 
empty(p_list list){
	while(list->next !=NULL){
		pop_elem_list(list);
	}
}

