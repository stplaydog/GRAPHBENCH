#include "graph.h"
#include "csr.h"

void read_graph_csr(t_csr *graph, char *file, char *dir, sint type){
	FILE *stream;
	char file_buf[100];
	char str[100];
	lint v_num;
	lint e_num;
	lint v_id, v_to;
	lint offset;
	lint e_weight=0;
	lint i,j;
	lint *v_idx_tmp;
	lint e_num_read = 0;
	lint g_v_size;
	
	sprintf(file_buf, "%s%s", dir, file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!", file);
		ERROR_PRINT();
	}
	graph->edge_idx = (lint*)_mm_malloc(sizeof(lint)*graph->e_size, 64);
	graph->edge_info = (t_e*)_mm_malloc(sizeof(t_e)*graph->e_size, 64);
	while(!feof(stream))
	{
		fscanf(stream, "%s", str);
		if(str[0]=='p' && str[1]=='\0'){
			fscanf(stream, "%s %lld %lld %lld %lld", str, &v_num, &e_num, &offset, &g_v_size);
			v_idx_tmp = (lint*) _mm_malloc(sizeof(lint)*v_num, 64);
			v_idx_tmp[0] = 0;	
			for(i=1;i<v_num;i++)
				v_idx_tmp[i] = graph->vet_idx[i-1];
		}
		else if(str[0]=='a' && str[1]=='\0'){
			fscanf(stream, "%lld %lld %lld",&v_id , &v_to, &e_weight);
			//graph->edge_info[e_num_read].edge_weight = (double)e_weight/(double)v_num;
			//if(v_to%2==0)
			//	graph->edge_info[e_num_read].edge_weight = 0.1;
			//else
			//	graph->edge_info[e_num_read].edge_weight = 0.2;
			graph->edge_info[e_num_read].edge_weight = (double)v_to/(double)v_num;
			if(type == RECV){
				v_id-=1;
				v_to-=offset;
			}
			else if(type == SEND){
				v_id-=offset;
				v_to-=1;
			}
			if(type==SEND){
				graph->edge_idx[v_idx_tmp[v_id]]=v_to;
				v_idx_tmp[v_id]++;
			}
			else if(type==RECV){
				graph->edge_idx[v_idx_tmp[v_to]]=v_id;
				v_idx_tmp[v_to]++;
			}
			e_num_read++;
			if(e_num_read == e_num) 
				break;	
		}
	}
	_mm_free(v_idx_tmp);
	fclose(stream);
}

void 
read_graph_csr_par(t_csr *graph, char *file, 
			char *dir, sint type)
{
}

void
separate_csr(t_csr *graph, t_csr *gl, t_csr *gh, double delta){
	gl->vet_idx = (lint*)_mm_malloc(sizeof(lint)*graph->v_size, 64);		
	gh->vet_idx = (lint*)_mm_malloc(sizeof(lint)*graph->v_size,64);		
	gl->vet_info = (t_vi*)_mm_malloc(sizeof(t_vi)*graph->v_size, 64);		
	gh->vet_info = (t_vi*)_mm_malloc(sizeof(t_vi)*graph->v_size,64);		
	gl->v_size = gh->v_size = graph->v_size;
	gl->e_size = gh->e_size = 0;
	int vet, idx=0;
	for(vet=0;vet<graph->v_size;vet++){
		for(; idx<graph->vet_idx[vet]; idx++){
			if(graph->edge_info[idx].edge_weight<=delta){
				gl->vet_info[vet].vet_deg++;
				gl->e_size++;
			}	
			else{
				gh->vet_info[vet].vet_deg++;
				gh->e_size++;
			}
		}
	}
	gl->edge_idx = (lint*)_mm_malloc(sizeof(lint)*gl->e_size, 64);
	gh->edge_idx = (lint*)_mm_malloc(sizeof(lint)*gh->e_size, 64);
	gl->edge_info = (t_e*)_mm_malloc(sizeof(t_e)*gl->e_size, 64);
	gh->edge_info = (t_e*)_mm_malloc(sizeof(t_e)*gh->e_size, 64);
	lint i;
	lint sum_l=0, sum_h=0;
	for(i=0;i<graph->v_size; i++){
		sum_l += gl->vet_info[i].vet_deg;
		gl->vet_idx[i]=sum_l;
		sum_h += gh->vet_info[i].vet_deg;
		gh->vet_idx[i]=sum_h;
	}	
	lint *v_idx_tmp_l = (lint*) _mm_malloc(sizeof(lint)*gl->v_size, 64);
	lint *v_idx_tmp_h = (lint*) _mm_malloc(sizeof(lint)*gh->v_size, 64);
	v_idx_tmp_l[0] = 0;	
	v_idx_tmp_h[0] = 0;	
	for(i=1;i<gl->v_size;i++)
		v_idx_tmp_l[i] = gl->vet_idx[i-1];
	for(i=1;i<gh->v_size;i++)
		v_idx_tmp_h[i] = gh->vet_idx[i-1];
	idx=0;
	for(vet=0;vet<graph->v_size;vet++){
		for(; idx<graph->vet_idx[vet]; idx++){
			if(graph->edge_info[idx].edge_weight<=delta){
				gl->edge_idx[v_idx_tmp_l[vet]] = graph->edge_idx[idx];	
				gl->edge_info[v_idx_tmp_l[vet]].edge_weight = graph->edge_info[idx].edge_weight;	
				v_idx_tmp_l[vet]++;
			}
			else{
				gh->edge_idx[v_idx_tmp_h[vet]] = graph->edge_idx[idx];	
				gh->edge_info[v_idx_tmp_h[vet]].edge_weight = graph->edge_info[idx].edge_weight;	
				v_idx_tmp_h[vet]++;
			}
		}
	}
}


void scan_csr_idx(t_csr *graph, char *file, char *dir, sint type)
{
	FILE *stream;
	char file_buf[100];
	char str[100];
	lint v_num;
	lint e_num;
	lint offset;
	lint g_v_size;
	lint count=0;

	lint v_id=0;
	lint v_to=0;
	lint e_weight=0;
	lint i,j,sum=0;
	lint e_num_read = 0;

	sprintf(file_buf, "%s%s", dir, file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		ERROR_PRINT();
	}
	while(!feof(stream))
	{
		fscanf(stream, "%s", str);
		if(str[0]=='p' && str[1]=='\0'){
			fscanf(stream, "%s %lld %lld %lld %lld\n", str, &v_num, &e_num, &offset, &g_v_size);
			//printf("v_num: %d e_num: %d\n", v_num, e_num);
			graph->offset = offset;
			graph->v_size = v_num;
			graph->e_size = e_num;
			graph->global_v_size = g_v_size;
			graph->vet_idx = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
			graph->vet_info = (t_vi*)_mm_malloc(sizeof(t_vi)*v_num, 64);
			for(i=0;i<v_num;i++){
				graph->vet_info[i].weight = 1.0;
				graph->vet_info[i].wb = 0;
				graph->vet_info[i].recip = 0;
				graph->vet_info[i].vet_deg = 0;
			}
		}
		else if(str[0]=='a' && str[1]=='\0'){
			fscanf(stream, "%lld %lld %lld\n",&v_id , &v_to, &e_weight);
			if(type == RECV){
				v_id-=1;
				v_to-=offset;
			}
			else if(type == SEND)
			{
				v_id-=offset;
				v_to-=1;
			}
			if(type==SEND)
				graph->vet_info[v_id].vet_deg++;		
			else if(type == RECV)
				graph->vet_info[v_to].vet_deg++;
			e_num_read++;
			if(e_num_read == e_num) 
				break;
		}
	}
	for(i=0;i<v_num;i++){
		sum += graph->vet_info[i].vet_deg;
		graph->vet_info[i].recip=1.0/(double)graph->vet_info[i].vet_deg;
		//if(i==0)
		//	printf("transfered deg is %1.10f, recip is %1.10f\n", (double)graph->vet_info[i].vet_deg, graph->vet_info[i].recip);
		graph->vet_idx[i]=sum;
	}
	fclose(stream);
}


	void 
scan_csr_idx_par(t_csr *graph, char *file, 
		char *dir, sint type)
{
}

void 
read_csr_bin(t_csr *graph, char* file,
		char* dir, sint type, int myrank){
	lint i;
	char f_buf[100];
	FILE *reader;
	sprintf(f_buf, "%s%s",dir, file);
	if((reader = fopen(f_buf, "rb"))==NULL){
		DPRINTF(1, "the file %s does not exist\n", f_buf);
		ERROR_PRINT();
	}
	lint v_num, e_num, offset, g_v_size;
	fread(&v_num, sizeof(lint), 1, reader);
	fread(&e_num, sizeof(lint), 1, reader);
	fread(&offset, sizeof(lint), 1, reader);
	fread(&g_v_size, sizeof(lint), 1, reader);
	lint *edge_arr = (lint*)_mm_malloc(sizeof(lint)*e_num*2, 64);
	fread(edge_arr, sizeof(lint), e_num*2, reader);
	fclose(reader);
	//init graph
	graph->offset = offset;
	graph->v_size = v_num;
	graph->e_size = e_num;
	graph->global_v_size = g_v_size;
	graph->vet_idx = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	graph->vet_info = (t_vi*)_mm_malloc(sizeof(t_vi)*v_num, 64);
	graph->edge_idx = (lint*)_mm_malloc(sizeof(lint)*e_num, 64);
	graph->edge_info = (t_e*)_mm_malloc(sizeof(t_e)*graph->e_size, 64);
	//MPI_Barrier(MPI_COMM_WORLD);
	for(i=0;i<v_num;i++){
		graph->vet_info[i].weight = INFINITY;
		graph->vet_info[i].wb = 0;
		graph->vet_info[i].recip = 0;
		graph->vet_info[i].vet_deg = 0;
	}
	//scan the edge list
	lint *tmp = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	for(i=0;i<v_num;i++)
		tmp[i]=0;

	for(i=0;i<e_num;i++){
		lint from = edge_arr[i*2];
		lint to = edge_arr[i*2+1];
		if(type==RECV){
			to-=offset;
			tmp[to]++;
		}
		else if(type == SEND){
			from -= offset;
			tmp[from]++;
		}
	}
	//printf("Node: %d done actual reading file %s, begin processing\n", myrank, f_buf);
	//MPI_Barrier(MPI_COMM_WORLD);
	lint sum=0;
	for(i=0;i<v_num;i++){
		graph->vet_idx[i]=tmp[i]+sum;
		sum+=tmp[i];
	}
	lint *tmp_idx = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	tmp_idx[0]=0;
	for(i=1;i<v_num;i++){
		tmp_idx[i]=graph->vet_idx[i-1];
		//if(myrank==1)
		//printf("i: %lld idx: %lld myrank: %d v_num:%lld\n", i, tmp_idx[i], myrank,v_num);
	}
	for(i=0;i<e_num;i++){
		lint from = edge_arr[i*2];
		lint to = edge_arr[i*2+1];
		if(type==RECV){
			to -= offset;
			from -= OFFSET;
			graph->edge_idx[tmp_idx[to]]=from;
			tmp_idx[to]+=1;
		}
		else if(type == SEND){
			from -= offset;
			to -= OFFSET;
			graph->edge_idx[tmp_idx[from]]=to;
			tmp_idx[from]+=1;
		}
	}
	for(i=0;i<e_num;i++){
		if(i%2==0)
			graph->edge_info[i].edge_weight = 0.1;	
		else
			graph->edge_info[i].edge_weight = 0.2;
	}
}

void vis_graph_csr(t_csr *graph)
{
	lint idx=0;
	lint vet;
	FILE *f;
	if((f = fopen("./vis/csr.vis", "w"))==NULL){
		printf("the file you input does not exist!");
		ERROR_PRINT();
	}

	fprintf(f, "digraph G {\n");
	for(vet=0;vet<graph->v_size;vet++){
		for(; idx<graph->vet_idx[vet]; idx++){
			fprintf(f, "%lld -> %lld;\n", vet, graph->edge_idx[idx]);
		}
	}
	fprintf(f, "}\n");
	fclose(f);
}

//void list_path(char *file, char *dir, lint target){
//	t_csr *gr = (t_csr*)malloc(sizeof(t_csr));
//	scan_csr_idx(gr, file, dir, RECV);
//	read_graph_csr(gr, file, dir, RECV);
//	p_list list = (p_list)malloc(sizeof(t_list));	
//	creat_list(list);
//	add_list_one(list, source);
//	
//}
