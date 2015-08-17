#include "graph.h"
#include "list.h"

void
print_B(p_list *B, FILE *log, lint iter);
void
relax(t_csr *gs, p_list *precedent,
		p_list *B, p_list R, 
		double delta, lint *b_total);
lint
find_edge_idx(t_csr *gs, lint source, lint target);
void
request(t_csr *gs, double delta, 
		p_list list, p_list R, 
		sint type);
void
remember(p_list B, p_list S, lint *b_total);
