#include "graph.h"
#include "csr_sspath.h"


void run_djkstra(char *file, char *dir, lint source){
	t_csr *gs = (t_csr*)malloc(sizeof(t_csr));
	read_csr_bin(gs, file,  dir, SEND, 0);
	lint *precedent = (lint*)malloc(sizeof(lint)*gs->v_size);
	int i;
	gs->vet_info[source].weight = 0.0;
	p_list Q = (p_list)malloc(sizeof(t_list));
	add_list_two(Q, source, 0.0);	
	lint iter=0;
	while(Q->next != NULL){
		p_list elem = pop_min(Q); 
		lint from = elem->int_val;
		double from_val = elem->double_val;
		lint start = from == 0 ? 0:gs->vet_idx[from-1];
		lint end = gs->vet_idx[from];
		for(i=start; i<end; i++){
			lint to = gs->edge_idx[i];
			double edge_weight = gs->edge_info[i].edge_weight;	
			double to_val = gs->vet_info[to].weight; 
			if((from_val+edge_weight)<to_val){
				//DPRINTF(1, "%f %f %f\n", from_val, edge_weight, to_val);
				gs->vet_info[to].weight = from_val+edge_weight;	
				precedent[to] = from;
				add_list_two(Q, to, gs->vet_info[to].weight);
			}
		}
		//DPRINTF(1, "this is iter %lld\n", iter++);
	}
	for(i=0;i<gs->v_size;i++)
		DPRINTF(1, "%f ", gs->vet_info[i].weight);
	DPRINTF(1, "\n");
}

void
run_sspath(t_csr *gs, p_list *predecesor, 
		lint source, double delta){
	int i;	
	char f_buf[100];	
	sprintf(f_buf, "./data/log/ori");
	FILE *log = fopen(f_buf, "w");	
	p_list *B = (p_list*)malloc(sizeof(p_list)*1000);		
	for(i=0;i<1000;i++){
		B[i] = (p_list) malloc(sizeof(t_list));
		creat_list(B[i]);
	}
	p_list R = (p_list)malloc(sizeof(t_list));
	p_list S = (p_list)malloc(sizeof(t_list));

	creat_list(R);
	creat_list(S);
	for(i=0;i<gs->v_size;i++)
		gs->vet_info[i].weight = INFINITY;
	//relax source vertex
	add_list_two(B[0], source, 0.0);
	gs->vet_info[source].weight = 0.0;
	lint b_total = 1;
	lint iter = 0;

	//major computing
	while(b_total>0){
		empty(R);	
		while (B[iter]->next!=NULL){
			request(gs, delta, B[iter], R, TYPE_LIGHT);
			remember(B[iter], S, &b_total);
			relax(gs, predecesor, B, R, delta, &b_total);
			//DPRINTF(1, "finished relax\n");
		}
		request(gs, delta, S, R, TYPE_HEAVY);
		relax(gs, predecesor, B, R, delta, &b_total);
		print_B(B, log, iter);
		//printf("=================iter %ld\n", iter);
		iter++;
		//DPRINTF(1, "finished iteration %lld \n", iter);
	}
	//DPRINTF(1, "iterations taken %lld\n", iter);	
	//for(i=0;i<gs->v_size;i++)
	//	DPRINTF(1, "vertex %lld %f\n", (i+1), gs->vet_info[i].weight);
	//DPRINTF(1, "\n");
}

void
relax(t_csr *gs, p_list *predecesor, 
		p_list *B, p_list R, 
		double delta, lint *b_total){
	p_list tmp = R->next;
	while(tmp!=NULL){
		lint source_idx = tmp->int_val;
		double source_val = gs->vet_info[source_idx].weight;
		lint target_idx = tmp->int_val1;
		lint edge_idx = find_edge_idx(gs, source_idx, target_idx);
		double target_ori_val = gs->vet_info[target_idx].weight;
		double target_val = gs->edge_info[edge_idx].edge_weight;
		double target_new_val = target_val+source_val;
		//if(target_idx == 48)
		//	printf("source %lld target %lld o_val %lf new_val %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
		//DPRINTF(1, "%f %f\n", target_ori_val, target_new_val);
		if(target_new_val <= target_ori_val){
			gs->vet_info[target_idx].weight = target_new_val;
			if(target_ori_val!=INFINITY){
				lint old_buck = (lint)(target_ori_val/delta);
				lint new_buck = (lint)(target_new_val/delta);	
				if(old_buck!=new_buck){
					pop_find(B[old_buck], target_idx);
					add_list_two(B[new_buck], target_idx, gs->vet_info[target_idx].weight);
				}
			}
			else{
				lint new_buck = (lint)(target_new_val/delta);
				add_list_two(B[new_buck], target_idx, gs->vet_info[target_idx].weight);
				(*b_total)++;
			}
			if(target_new_val < target_ori_val)
				empty(predecesor[target_idx]);	
			add_list_one(predecesor[target_idx], source_idx);
		}
		tmp = tmp->next;
	}	
}

lint
find_edge_idx(t_csr *gs, lint source, lint target){
	lint start = source==0?0:gs->vet_idx[source-1];
	lint end = gs->vet_idx[source];
	lint i;
	for(i=start; i<end; i++){
		if(gs->edge_idx[i]==target)
			return i;
	}
	return -1;
}

void
request(t_csr *gs, double delta, 
		p_list list, p_list R, 
		sint type){
	p_list tmp = list->next;
	while(tmp!=NULL){
		lint query = tmp->int_val;
		lint start = query==0?0:gs->vet_idx[query-1];	
		lint end = gs->vet_idx[query];
		lint i;
		for(i=start; i<end; i++){
			lint target = gs->edge_idx[i];
			double target_val = gs->edge_info[i].edge_weight;
			if((target_val<=delta && type == TYPE_LIGHT) ||
					(target_val>delta && type == TYPE_HEAVY)){
				add_list_three(R, query, target, target_val);
			}
		}
		tmp = tmp->next;
	}
}

void
remember(p_list B, p_list S, lint *b_total){
	p_list tmp = pop_elem_list(B);
	while(tmp!= NULL){
		lint vet_id = tmp->int_val;
		add_list_one(S, vet_id);
		tmp = pop_elem_list(B);
		(*b_total)--;
	}
}

void
print_B(p_list *B, FILE *log, lint iter){
	lint i, j, k;
	fprintf(log, "=============After %lld, the B Map is:==================\n", iter);
	for(i=0;i<buck_size; i++){
		if(B[i]->next!= NULL)
			fprintf(log, "Map on iter %lld\n", i);
		else 
			continue;
		p_list tmp = B[i]->next;
		while(tmp != NULL){
			//printf("%lld|%lld|%lf ", tmp->int_val, tmp->int_val1, tmp->double_val);
			fprintf(log, "%lld ", tmp->int_val);
			tmp = tmp->next;
		}
		fprintf(log, "\n");
	}
}
