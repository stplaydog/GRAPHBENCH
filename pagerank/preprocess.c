#include "preprocess.h"

void
reduce_graph_for_pr(char *file, char *dir){
	lint i,j,k;
	sint myrank, size;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	
	char f_buf[100];
	sprintf(f_buf, "%sdatex/%s_sd_%d_%d", dir, file, myrank, num_nodes);
	FILE *reader;
	if((reader=fopen(f_buf, "r"))==NULL){
		printf("file %s not exist\n", f_buf);
	}

	//read graph file--------------------------	
	lint num_v, num_e, offset, g_v_num;
	fscanf(reader, "p sp %lld %lld %lld %lld\n", &num_v, &num_e, &offset, &g_v_num);
	lint *arr_from = (lint*)_mm_malloc(sizeof(lint)*num_e, 64);
	lint *arr_to = (lint*)_mm_malloc(sizeof(lint)*num_e, 64);
	lint num_e_read=0;
	lint from, to, e_weight;
	while(1){
		fscanf(reader, "a %lld %lld %lld\n", &from, &to, &e_weight);
		arr_from[num_e_read] = from;
		arr_to[num_e_read] = to -1;
		num_e_read++;
		if(num_e_read==num_e)
			break;
	}
	fclose(reader);

	//read send ub-----------------------------
	sprintf(f_buf, "%sparex/%s_send_%d", dir, file, num_nodes);
	FILE *reader_sdub;
	if((reader_sdub=fopen(f_buf, "r"))==NULL){
		printf("file %s not exist\n", f_buf);
	}
	lint **send_ub = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	//send_ub_new is used to get the reduced send ub
	lint **send_ub_new = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	for(i=0;i<num_nodes;i++){
		send_ub[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		send_ub_new[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		for(j=0;j<num_nodes;j++)
			fscanf(reader_sdub, "%lld", &(send_ub[i][j]));
	}
	fclose(reader_sdub);
	
	//read parition----------------------------
	sprintf(f_buf, "%sparex/%s_recip_%d_%d", dir, file, myrank, num_nodes);
	FILE *reader_recip;
	if((reader_recip=fopen(f_buf, "r"))==NULL){
		printf("file %s not exist\n", f_buf);
	}
	lint num_v_read=0;
	double *recip = (double*)_mm_malloc(sizeof(double)*num_v, 64);
	while(1){
		fscanf(reader_recip, "%lf", &recip[num_v_read]);
		num_v_read++;
		if(num_v_read == num_v)
			break;
	}
	lint *partition = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
	for(i=0;i<num_nodes;i++){
		fscanf(reader_recip, "%lld", &partition[i]);
	}
	fclose(reader_recip);
	
	//assign value to the arrays and sort them---------------------
	lint **parr = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64); //used to store the real value ex: 2 -> 3 store 3
	lint **parr_idx = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64); //used to store the index of the value in the parr
	lint **parr_idx_map = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64); //
	lint **parr_rev_idx = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	lint *par_start = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
	lint *idx = (lint*)_mm_malloc(sizeof(lint)*num_e, 64); //idx -> parr_idx and know the final position for a given vertex of a given edge
	for(i=0;i<num_nodes;i++){
		parr[i] = (lint*)_mm_malloc(sizeof(lint)*send_ub[myrank][i], 64);
		parr_idx[i] = (lint*)_mm_malloc(sizeof(lint)*send_ub[myrank][i], 64);
		parr_idx_map[i] = (lint*)_mm_malloc(sizeof(lint)*send_ub[myrank][i], 64);
		parr_rev_idx[i] = (lint*)_mm_malloc(sizeof(lint)*send_ub[myrank][i], 64);
		par_start[i] = 0;
	}
	for(i=0; i<num_e; i++){
		lint node_id = bin_search(partition, arr_to[i], 0, num_nodes-1); 	
		lint pos = par_start[node_id];
		parr[node_id][pos] = arr_to[i];
		parr_idx[node_id][pos] = pos;
		idx[i] = par_start[node_id];
		par_start[node_id]++;
	}
	for(i=0;i<num_nodes;i++)
		q_sort_two(parr[i], parr_idx[i], 0, send_ub[myrank][i]-1);
	MPI_Barrier(MPI_COMM_WORLD);
	
	//reduce the redundant destination vertices-------------------
	for(i=0;i<num_nodes;i++){
		lint sum=0;
		parr_idx_map[i][0] = 0;
		for(j=1;j<send_ub[myrank][i];j++){
			if(parr[i][j] == parr[i][j-1])
				parr_idx_map[i][j] = sum;
			else{
				parr_idx_map[i][j] = ++sum;
			}
			parr_rev_idx[i][parr_idx[i][j]]= parr_idx_map[i][j];
		}
		send_ub_new[myrank][i] = sum;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//write back-------------------------------------------------
	sprintf(f_buf, "%sparex/%s_pidx_%d_%d", dir, file, myrank, num_nodes);
	FILE *writer_pidx;
        if((writer_pidx=fopen(f_buf, "w"))==NULL){
                printf("file %s not exist\n", f_buf);
        }
	for(i=0; i<num_nodes; i++)
		par_start[i]=0;
	for(i=0;i<num_e;i++){
		lint node_id = bin_search(partition, arr_to[i], 0, num_nodes-1);
		lint pos = par_start[node_id];
		lint index = parr_rev_idx[node_id][pos];
		par_start[node_id]++;
		fprintf(writer_pidx, "%lld\n", index);
	}
	fprintf(writer_pidx, "\n");
	fclose(writer_pidx);
	
	//write new send ub-------------------------------------------
	if(myrank == 0){
		for(i=1;i<num_nodes;i++){
			MPI_Recv(send_ub_new[i], num_nodes, MPI_LONG, i, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	else
		MPI_Send(send_ub_new[myrank], num_nodes, MPI_LONG,  0, 1, MPI_COMM_WORLD );
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0){
		sprintf(f_buf, "%sparex/%s_send_%d", dir, file, num_nodes);
		FILE *writer_sdub;
		if((writer_sdub=fopen(f_buf, "w"))==NULL){
			printf("file %s not exist\n", f_buf);
		}
		for(i=0; i<num_nodes; i++){
			for(j=0;j<num_nodes; j++){
				fprintf(writer_sdub, "%lld ", send_ub_new[i][j]);
			}
			fprintf(writer_sdub, "\n");
		}
		fclose(writer_sdub);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
