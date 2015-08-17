#include "graph.h"
#include "csr_bc.h"
#include "list.h"
void
comp_num_spath(lint source, lint *num_spath, 
		lint *dis_idx, p_list *predecesor, 
		lint v_size);
void 
comp_delta(lint source, lint *num_spath, 
		double *delta, lint *dis_idx, 
		p_list *predecesor, lint v_size);
void 
empty_predecesor(p_list *predecesor, lint v_size);


	void 
run_bc(char *file, char *dir)
{
	lint i, j;
	t_csr *gs = (t_csr*)malloc(sizeof(t_csr)); 	
	if(bin==FALSE){
		scan_csr_idx(gs, file, dir, SEND);
		read_graph_csr(gs, file, dir, SEND);
		//vis_graph_csr(gs);
	}
	else
		read_csr_bin(gs, file, dir, SEND, 0);
	p_list *predecesor = (p_list *)malloc(sizeof(p_list)*gs->v_size);
	double *delta = (double*)malloc(sizeof(double)*gs->v_size);	
	double *total_delta = (double*)malloc(sizeof(double)*gs->v_size);	
	double *distance = (double*)malloc(sizeof(double)*gs->v_size);
	lint *dis_idx = (lint*)malloc(sizeof(lint)*gs->v_size);
	lint *num_spath = (lint*)malloc(sizeof(lint)*gs->v_size);
	for(i=0;i<gs->v_size;i++){
		predecesor[i] = (p_list)malloc(sizeof(t_list));
		creat_list(predecesor[i]);
		total_delta[i]=0;
	}
	// brandes algorithm
	for(i=0;i<gs->v_size;i++){
		//printf("processing vertex %lld\n", i);
		int source = i;
		run_sspath(gs, predecesor, source, 0.1);
		for(j=0;j<gs->v_size;j++){
			distance[j] = gs->vet_info[j].weight;
			dis_idx[j] = j;
		}
		q_sort_double(distance, dis_idx, 0, gs->v_size-1);
	//	for(j=0;j<gs->v_size;j++)
	//		printf("%lf|%lld ", distance[j], dis_idx[j]);
	//	printf("\n");
		comp_num_spath(source, num_spath, dis_idx, predecesor, gs->v_size);
		//for(j=0;j<gs->v_size;j++)
		//	printf("%lld ", num_spath[j]);
		//printf("\n");
		comp_delta(source, num_spath, delta, dis_idx, predecesor, gs->v_size);
		for(j=0;j<gs->v_size;j++){
			if(j==source)
				continue;
			total_delta[j] += delta[j];
		}
		empty_predecesor(predecesor, gs->v_size);
	}
	printf("the final result for betweeness centrality is:\n");
	for(i=0;i<gs->v_size; i++)
		printf("%lld:%lf ", (i+1), total_delta[i]);
	printf("\n");
}

void
comp_num_spath(lint source, lint *num_spath, 
		lint *dis_idx, p_list *predecesor, 
		lint v_size){
	lint i;
	for(i=0; i<v_size; i++){
		num_spath[i]=0;
	}
	num_spath[source]=1;
	for(i=0; i<v_size; i++){
		if(dis_idx[i]==source)
			continue;
		lint w = dis_idx[i];		
		//printf("source %lld w %lld\n",source, w);
		p_list tmp = predecesor[w]->next;
		while(tmp != NULL){
			lint v = tmp->int_val; 
			//if(source==0)
			//	printf("w %lld val: %lld v %lld val %lld\n", w, num_spath[w], v, num_spath[v]);
			num_spath[w]+=num_spath[v];
			tmp = tmp->next;
		}
	}
//	for(i=0;i<v_size;i++)
//		printf("%lld:%lld ", (i+1), num_spath[i]);
//	printf("\n");
}

void 
comp_delta(lint source, lint *num_spath, 
		double *delta, lint *dis_idx, 
		p_list *predecesor, lint v_size){
	lint i;
	for(i=0;i<v_size;i++){
		delta[i]=0.0;
	}
	for(i=v_size-1;i>=0;i--){
		if(dis_idx[i]==source)
			continue;
		lint w = dis_idx[i];	
		p_list tmp = predecesor[w]->next;
		while(tmp != NULL){
			lint v = tmp->int_val;
			//if(source==15)
			//	printf("v: %lld w: %lld v_num %lld w_num %lld\n", v, w, num_spath[v], num_spath[w]);
			if(v!=source)
				delta[v]+=(1.0+delta[w])*(double)num_spath[v]/(double)num_spath[w];
			tmp = tmp->next;
		}	
	}
	//for(i=0;i<v_size;i++)
	//	printf("%lld:%lf ", (i+1), delta[i]);
	//printf("\n");
}

void 
empty_predecesor(p_list *predecesor, lint v_size){
	int i;
	for(i=0;i<v_size;i++)
		empty(predecesor[i]);	
}
