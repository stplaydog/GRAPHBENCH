#include "graph.h"
#include "csr_partition.h"
#include "utils.h"

void comp_distrib_csr(t_csr *gs, t_csr *gr, 
			lint *distrib, lint *send_size, 
			lint *recv_size, sint type)
{
	lint distrib_num = gr->e_size/num_nodes;
	lint i,j=0;
	for(i=0;i<num_nodes;i++)
	{
		distrib[i]=0;
		send_size[i]=0;
		recv_size[i]=0;
	}
	for (i=0;i<gs->v_size;i++)
	{
		send_size[j] += gs->vet_info[i].vet_deg;
		recv_size[j] += gr->vet_info[i].vet_deg;
		DPRINTF(3, "%lld|%lld ", i, gs->vet_info[i].vet_deg);
		if((send_size[j]>distrib_num && type==SEND) 
			|| (recv_size[j]>distrib_num && type == RECV))
		{
			distrib[j]=i;
			j++;
		}
	}
	distrib[num_nodes-1]=gr->v_size-1;
}

void chunk_files(char* file, char *dir, 
			lint *distrib, lint *send_size, 
			lint *recv_size)
{
	lint i,j=0;
	char name_buf_send[100];
	char name_buf_recv[100];
	FILE *f_in;
	FILE *writer_send[num_nodes];
	FILE *writer_recv[num_nodes];
	char str[100];
	lint v_num, e_num, offset;
	lint v_idx, v_to, e_weight;
	char extension[100];
	lint num_read =0;
	
	lint **send_ub = (lint**)malloc(sizeof(lint*)*num_nodes);
	for(i=0;i<num_nodes;i++)
	{
		send_ub[i] = (lint*)malloc(sizeof(lint)*num_nodes);
		for(j=0; j<num_nodes; j++)
			send_ub[i][j]=0;
	}
	FILE *send_f;
	char send_buf[100];
	sprintf(send_buf, "%s%s%s_send_%s", dir, par_extend, file, nodes_extension);
	if((send_f = fopen(send_buf, "w"))==NULL){
		printf("the file %s does not exists!\n", send_buf);
		ERROR_PRINT();
	}

	for(i=0;i<num_nodes;i++){
		snprintf(extension, 10,"%lld",i);
		concate(name_buf_send, 7, dir, data_extend, file, "_sd_", extension, "_", nodes_extension);
		concate(name_buf_recv, 7, dir, data_extend, file, "_rc_", extension, "_", nodes_extension);

		if((writer_send[i]=fopen(name_buf_send,"w"))==NULL){
			printf("the file you input does not exist!");
			ERROR_PRINT();
		}
		if((writer_recv[i]=fopen(name_buf_recv,"w"))==NULL){
			printf("the file you input does not exist!");
			ERROR_PRINT();
		}
		DPRINTF(3, "open file %s\n", name_buf_send);

	}
	char f_buf[100];
	sprintf(f_buf, "%s%s",dir, file);
	if((f_in=fopen(f_buf,"r"))==NULL){
		printf("the file you input does not exist!");
		ERROR_PRINT();
	}
	while(!feof(f_in)){
		fscanf(f_in, "%s", str);
		if(str[0]=='p' && str[1]=='\0'){
			fscanf(f_in, "%s %lld %lld %lld", str, &v_num, &e_num, &offset);
			for(i=0;i<num_nodes; i++){
				lint size = (i==0?(distrib[0]+1):distrib[i]-distrib[i-1]);
				lint off = (i==0?0:distrib[i-1]+1)+1;
				fprintf(writer_send[i], "p sp %lld %lld %lld %lld\n", size, send_size[i], off, v_num);
				fprintf(writer_recv[i], "p sp %lld %lld %lld %lld\n", size, recv_size[i], off, v_num);
			}
		}
		if(str[0]=='a' && str[1]=='\0'){
			fscanf(f_in, "%lld %lld %lld", &v_idx, &v_to, &e_weight);
			v_idx -= 1;
			v_to -= offset;
			lint send_id = bin_search(distrib, v_idx, 0, num_nodes-1);
			lint recv_id = bin_search(distrib, v_to, 0, num_nodes-1);
			if(send_id!=recv_id)
				send_ub[send_id][recv_id]++;	
			fprintf(writer_send[send_id], "a %lld %lld %lld\n", v_idx+1, v_to+offset, e_weight);	
			fprintf(writer_recv[recv_id], "a %lld %lld %lld\n", v_idx+1, v_to+offset, e_weight);	
			num_read++;
			if(num_read==e_num)
				break;
		}
	}
	for(i=0;i<num_nodes;i++)
	{
		for(j=0;j<num_nodes;j++)
		{
			fprintf(send_f, "%lld ", send_ub[i][j]);
		}
		fprintf(send_f, "\n");
	}
	for(i=0;i<num_nodes;i++){
		fclose(writer_send[i]);
		fclose(writer_recv[i]);
	}
	fclose(f_in);
	fclose(send_f);
}

void init_node_id(t_csr *graph, lint *distrib, sint num_nodes){
	lint i;
	for(i=0;i<graph->e_size;i++)
		graph->edge_info[i].node_id = bin_search(distrib, graph->edge_idx[i], 0, num_nodes-1); 	
}

void write_recip(char* file, char* dir, 
		lint *distrib, t_csr *graph_send)
{
	lint i,j=0,k;
	FILE *writer[num_nodes];
	char f_buf[100];
	char extension[100];
	for(i=0;i<num_nodes;i++){
		sprintf(f_buf, "%s%s%s_recip_%lld_%s", dir, par_extend, file, i, nodes_extension);
		if((writer[i] = fopen(f_buf, "w"))==NULL){
			printf("the file does not exists!\n");
			ERROR_PRINT();
		}
	}
	for(i=0;i<graph_send->v_size;i++){
		fprintf(writer[j], "%f ", graph_send->vet_info[i].recip);
		if(i<distrib[j])
			;
		else if(i==distrib[j])
			j++;
	}
	for(i=0;i<num_nodes;i++){
		fprintf(writer[i], "\n");
		for(k=0;k<num_nodes;k++)
			fprintf(writer[i], "%lld ", distrib[k]);
		fprintf(writer[i], "\n");
	}
	for(i=0;i<num_nodes;i++){
		fprintf(writer[i], "\n");
		fclose(writer[i]);
	}
}

void read_recip(char* file, char *dir, lint *distrib, 
		t_csr *graph, sint myrank){
	lint i,j,k;
	FILE *reader;
	char f_buf[100];
	char extension[100];
	double in;
	sprintf(f_buf, "%s%s%s_recip_%d_%s", dir, par_extend, file, myrank, nodes_extension);
	sprintf(extension, "%s%s%s_recip_%d_%s_vis", dir, par_extend, file, myrank, nodes_extension);
	FILE *writer = fopen(extension, "w");
	if((reader = fopen(f_buf, "r"))==NULL){
		ERROR_PRINT();
	}
	for(i=0;i<graph->v_size;i++){
		fscanf(reader, "%lf", &graph->vet_info[i].recip);
	}
	for(k=0;k<num_nodes;k++)
		fscanf(reader, "%lld", &distrib[k]);
	fclose(reader);
}

void read_send_ub(char *file, char *dir, lint **send_ub)
{
	lint i,j;
	FILE *f;
	char buf[100];
	sprintf(buf, "%s%s%s_send_%s", dir, par_extend, file, nodes_extension);
	if((f=fopen(buf, "r"))==NULL){
		printf("the file %s does not exists!\n", file);
		ERROR_PRINT();
	}
	for(i=0;i<num_nodes;i++)
		for(j=0;j<num_nodes;j++)
			fscanf(f, "%lld", &(send_ub[i][j]));
	fclose(f);
}
