#include "graph.h"

struct vet_pos
{
	sint buck_id;
	lint buck_pos;
};

typedef struct vet_pos *p_vidx;
typedef struct vet_pos t_vidx;

struct arr_elem
{
	lint int_val;
	lint int_val1;
	double double_val;
};
typedef struct arr_elem t_arr_elem;
typedef struct arr_elem *p_arr_elem;

struct list_arr
{
	t_arr_elem *arr;
	lint size;
	lint idx;
	lint num;
};
typedef struct list_arr t_arr;
typedef struct list_arr *p_arr;

void
create_arr(p_arr arr, lint size);
void 
add_arr_one(p_arr arr, lint val1);
void 
add_arr_two(p_arr arr, lint val1, double dval1);
void 
add_arr_three(p_arr arr, lint val1, lint val2, double dval1);
void
pop_elem_arr(p_arr arr, p_arr_elem elem);
void 
remove_find_arr(p_arr arr, lint idx);
void 
empty_arr(p_arr arr);
void 
merge_sort(p_arr arr, p_arr tmp);
void 
partition(p_arr arr, p_arr tmp, lint * bin_idx);
