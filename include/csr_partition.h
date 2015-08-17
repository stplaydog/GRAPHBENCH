#include "graph.h"

//csr_partition.c
void comp_distrib_csr(t_csr *graph_send, t_csr *graph_recv, 
			lint *distrib, lint *send_size, 
			lint *recv_size, sint type);
void chunk_files(char* file, char *dir, 
			lint *distrib, lint *send_size, 
			lint *recv_size);
void write_recip(char* file, char* dir, 
			lint *distrib, t_csr *graph_send);
void read_recip(char* file, char *dir, lint *distrib, 
			t_csr *graph, sint myrank);
void read_send_ub(char *file, char *dir, lint **send_ub);
void init_node_id(t_csr *graph, lint *distrib, sint num_nodes);
