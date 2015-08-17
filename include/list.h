#include "graph.h"

struct list_elem
{
	lint int_val;
	lint int_val1;
	double double_val;
	struct list_elem *next;
};
typedef struct list_elem t_list;
typedef struct list_elem *p_list;

void
creat_stack(p_list stack);
void
add_list_one(p_list stack, lint val_int);
void
add_list_two(p_list stack, lint val_int, double val_double);
void
add_list_three(p_list stack, lint val_int, lint val1_int, double val_double);
p_list
pop_elem_list(p_list list);
p_list
pop_min(p_list list);
p_list
pop_find(p_list list, lint query);
void 
empty(p_list list);
sint 
find(p_list list, lint query);
void
run_sspath(t_csr *gs, p_list *predecesor, 
		lint source, double delta);
