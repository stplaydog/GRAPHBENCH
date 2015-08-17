#include "graph.h"
//csr_partition_mpi.c
void 
chunk_files_physical(char* file, char *dir, lint *global_info, sint myrank);
	void 
chunk_files_physical_to_bin(char* file, char *dir, lint *global_info, sint myrank);
void 
comp_distrib_csr_par(char *file, char *dir, 
			lint *distrib, lint *send_size, 
			lint *recv_size, lint *global_info, 
			sint type, sint myrank);
void 
chunk_files_par(char* file, char *dir, 
			lint *distrib, lint *send_size, 
			lint *recv_size, sint type, 
			lint *global_info, sint myrank);
void 
remove_files(char* file, char* dir, sint myrank);
void 
write_recip_par(char* file, char* dir, 
			lint *distrib, t_csr *graph_send);
void 
read_recip_par(char* file, char *dir, lint *distrib, 
		t_csr *graph, sint myrank);
