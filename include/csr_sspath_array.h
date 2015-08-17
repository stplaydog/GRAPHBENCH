#include "graph.h"
#include "array.h"
#include "list.h"
#include "timer.h"

void reorder_arr(p_arr arr, t_csr *gs, double *tmp);
void			
relax_arr_by_partition(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
relax_arr_single(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
request_arr_single(t_csr *gs, p_arr arr, p_arr R);
void
remember_arr_single(p_arr B, p_arr S, 
		lint *b_total, p_vidx vidx);

//this is for parallel code
void
print_sssp_results(p_time_list pl);
void 
compute_request_partition(p_arr R);
void 
print_B_map(FILE *log, p_arr *B, 
		lint iter);
void
set_vidx(p_arr B, p_vidx vidx, sint iter);
int
check_B_total_on_iter(p_arr **B, lint iter);
void 
compute_request_idx(sint thread_id, t_csr *gs, t_arr *arr, FILE *log);
void 
reduce_request_idx();
void
balance_request(p_arr RR, p_arr R, 
		sint thread_id);
lint
find_edge_idx_arr(t_csr *gs, lint source, lint target);
void
relax_arr(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
request_arr(t_csr *gs, p_arr arr, p_arr R);
void
remember_arr(p_arr B, p_arr S, lint *b_total, p_vidx vidx);
lint 
check_total(lint *b_total);
void recover_arr(p_arr *B, p_vidx vidx, 
		lint iter);
void
compute_RR_partition(p_arr RR);
