#include "graph.h"

//csr.c
void read_graph_csr(t_csr *graph, char *file, char *dir, sint type);
void scan_csr_idx(t_csr *graph, char *file, char *dir, sint type);
void print_pr_csr(t_csr *graph);
void vis_graph_csr(t_csr *graph);
void 
scan_csr_idx_par(t_csr *graph, char *file, 
			char *dir, sint type);
void 
read_graph_csr_par(t_csr *graph, char *file, 
			char *dir, sint type);
void 
read_csr_bin(t_csr *graph, char* file,
		char* dir, sint type, int myrank);
