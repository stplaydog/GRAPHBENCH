#include "graph.h"
//csr_partition_bin_mpi.c
void
remove_duplicate(char* file, char*dir, sint myrank);
void comp_distrib_csr_par_bin(char* file, char* dir, 
				lint *distrib, lint *send_size,
				lint *recv_size, sint type, sint myrank);
	void 
chunk_files_bin_par(char* file, char *dir, 
		lint *distrib, sint type, 
		sint myrank);
void 
remove_files(char* file, char* dir, sint myrank);
int
bin_to_msg(lint *edge_arr, sint *node_mark, 
		lint *msg, lint len, 
		sint node_id, lint start, lint size);
int
bin_from_msg(lint *edge_recv_arr, lint *msg, 
		lint len, lint start);
void 
bin_check_msg_size(lint *msg_size, sint *node_mark, 
			lint len);
void 
read_recip_bin(char* file, char *dir, lint *distrib, 
		t_csr *graph, sint myrank);
