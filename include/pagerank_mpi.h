#include "graph.h"

//csr_mpi.c
lint 
cpy_from_buff_to_msg(t_msg *msg, lint ***int_buffer, 
			double ***double_buffer, sint myrk, 
			lint to, sint idx_buff, 
			lint idx_msg);
void init_msg(t_msg *msg, lint **partition, sint myrank);
void to_msg(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank, 
		t_buff *buf);
void print_msg(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank);
void print_msg_file_send(FILE *f, double *diff, t_msg *msg, 
		t_csr *gs, sint myrank);
void print_msg_file_recv(FILE *f, double *diff, t_msg *msg, 
		t_csr *gs, sint myrank);
void print_sendub_file(FILE *f, lint**send_ub);
void from_msg(t_msg *msg, t_csr *gr, 
		sint myrank, lint **send_ub, 
		double purpose_jump, lint **idx_start,
		lint **idx_end);
void run_pagerank_csr_mpi(char *file, char *dir);
void 
compute_bucket_role(sint *bucket, sint *role, sint ii);
void
init_t_buff(t_buff *buf, int size);
//reduced message, there is no need to design a new from_msg which is compatible with the new to_msg_reduce

/*assign diff_pr to diff_tmp
 * */
void assign_diff(double *diff_tmp, lint *diff_assign,
			double *diff_pr, t_csr *gs, sint myrank);

/* to_msg_reduce is to assign msg to different destination msg buffer
 * */
void to_msg_reduce(double *diff_tmp, t_msg *msg,
			lint *pidx_tmp, lint *idx_tmp, 
			lint *bin_idx, lint *node_idx_tmp);
void to_msg_reduce_serial(double *diff_tmp, t_msg *msg,
			lint *pidx_tmp, lint *idx_tmp,
			lint *node_idx_tmpi, lint size);
void set_msg(t_msg *msg, lint **send_ub, sint myrank);

/*partition message is used to partition different destination vets to different bins
 * to avoid contension in the to_msg step
 * */
void partition_msg(lint *diff_assign, lint *pidx_tmp, 
			lint *idx_tmp, lint *bin_idx, 
			lint *node_idx_tmp, lint size);
void reverse_diff_assign(lint *diff_assign, lint e_size);
