#include "graph.h"
#include "csr_partition_mpi.h"
#include "utils.h"

/**
 * @brief chunk a graph file into multiple small graph files
 *
 * chunk files into num_nodes files
 * each file will have the same number of edges
 * 
 * @param file the name of the large graph file
 * @param dir the directory of the large graph file
 */
	void 
chunk_files_physical(char* file, char *dir, lint *global_info, sint myrank)
{
	if(myrank==0){
		lint i,j;
		char str[100];
		lint v_num, e_num, offset, g_v_size;
		lint v_id, v_to, e_weight;
		lint avg_e_size;
		lint num_read=0;
		lint now_write=0;
		lint ind_read[num_nodes];
		FILE *reader;
		FILE *writer[num_nodes];
		char r_buf[100];
		char w_buf[100];
		//reader file
		sprintf(r_buf, "%s%s", dir, file);
		if((reader = fopen(r_buf, "r"))==NULL){
			DPRINTF(1, "the file %s does not exist!\n", r_buf);
			ERROR_PRINT();	
		}
		//writer file
		for(i=0; i<num_nodes; i++){
			ind_read[i]=0;
			sprintf(w_buf, "%s%s%s_phy_%lld_%d", dir, data_extend, file, i, num_nodes);
			if((writer[i] = fopen(w_buf, "w"))==NULL){
				DPRINTF(1, "the file %s does not exist!\n", w_buf);
				ERROR_PRINT();	
			}
		}
		//do chunking work
		while(!feof(reader)){
			fscanf(reader, "%s", str);	
			if(str[0]=='p' && str[1]=='\0'){
				fscanf(reader, "%s %lld %lld %lld %lld", str, &v_num, &e_num, &offset, &g_v_size);
				avg_e_size = e_num/num_nodes;
			}
			else if(str[0]=='a' && str[1]=='\0'){
				fscanf(reader, "%lld %lld %lld",&v_id , &v_to, &e_weight);
				fprintf(writer[now_write], "a %lld %lld %lld\n", v_id , v_to, e_weight);	
				ind_read[now_write]++;
				if(ind_read[now_write]==avg_e_size && now_write!=(num_nodes-1))
				{
					now_write++;
				}
				num_read++;
				if(num_read==e_num){
					break;
				}
			}
		}
		fclose(reader);
		for(i=0;i<num_nodes;i++){
			global_info[i]=ind_read[i];
			fclose(writer[i]);
		}
		global_info[num_nodes]=v_num;
		global_info[num_nodes+1]=e_num;
	}
}

	void 
chunk_files_physical_to_bin(char* file, char *dir, lint *global_info, sint myrank)
{
	if(myrank==0){
		lint i,j;
		char str[100];
		lint v_num, e_num, offset, g_v_size;
		lint v_id, v_to, e_weight;
		lint avg_e_size;
		lint num_read=0;
		lint now_write=0;
		lint ind_read[num_nodes];
		FILE *reader;
		FILE *writer[num_nodes];
		char r_buf[100];
		char w_buf[100];
		//reader file
		sprintf(r_buf, "%s%s", dir, file);
		if((reader = fopen(r_buf, "r"))==NULL){
			DPRINTF(1, "the file %s does not exist!\n", r_buf);
			ERROR_PRINT();	
		}
		//writer file
		for(i=0; i<num_nodes; i++){
			ind_read[i]=0;
			sprintf(w_buf, "%s%s_%lld", dir, file, i);
			if((writer[i] = fopen(w_buf, "wb"))==NULL){
				DPRINTF(1, "the file %s does not exist!\n", w_buf);
				ERROR_PRINT();	
			}
		}
		//do chunking work
		lint **edge_list = (lint**)malloc(sizeof(lint*)*num_nodes);
		while(!feof(reader)){
			fscanf(reader, "%s", str);	
			if(str[0]=='p' && str[1]=='\0'){
				fscanf(reader, "%s %lld %lld %lld %lld", str, &v_num, &e_num, &offset, &g_v_size);
				avg_e_size = e_num/num_nodes;
				for(i=0;i<num_nodes;i++)
					edge_list[i]=(lint*)malloc(sizeof(lint)*(avg_e_size*2+10000));
			}
			else if(str[0]=='a' && str[1]=='\0'){
				fscanf(reader, "%lld %lld %lld",&v_id , &v_to, &e_weight);
				edge_list[now_write][ind_read[now_write]*2]=v_id-1;
				edge_list[now_write][ind_read[now_write]*2+1]=v_to-1;
				ind_read[now_write]++;
				if(ind_read[now_write]==avg_e_size && now_write!=(num_nodes-1))
				{
					now_write++;
				}
				num_read++;
				if(num_read==e_num){
					break;
				}
			}
		}
		for(i=0;i<num_nodes;i++){
			//DPRINTF(1, "this is file:%lld \n", i);
			//for(j=0;j<ind_read[i];j++)
			//	DPRINTF(1, "%lld %lld\n", edge_list[i][2*j], edge_list[i][2*j+1]);
			fwrite(&ind_read[i], sizeof(lint), 1, writer[i]);
			DPRINTF(1, "num egdes for rank %lld: %lld %lld\n", i, ind_read[i], e_num);
			fwrite(edge_list[i], sizeof(lint), ind_read[i]*2, writer[i]);
		}
		fclose(reader);
		for(i=0;i<num_nodes;i++){
			global_info[i]=ind_read[i];
			fclose(writer[i]);
		}
		global_info[num_nodes]=v_num;
		global_info[num_nodes+1]=e_num;
	}
}
/**
 * @brief read the structure of graph and partition it into subgraphs 
 *
 * compute the distribution of each node which has differnet number of 
 * vertices and euqal number of egdes
 *
 * @param file base file name
 * @param dir base directory
 * @param distrib how vertices are distribute among different files
 * @param send_size how many outgoing edges a subset of vertices has
 * @param recv_size how many ingoing edges a subset of vertices has
 * @param partition by #outgoing edges =>SEND #ingoing edges =>RECV
 */
	void 
comp_distrib_csr_par(char* file, char *dir, 
		lint *distrib, lint *send_size, 
		lint *recv_size, lint *global_info, 
		sint type, sint myrank)
{
	lint i=0,j=0;
	lint *v_num_send;
	lint *v_num_recv;
	lint *vg_num_send;
	lint *vg_num_recv;
	FILE *reader;	
	char f_buf[100];
	char str[100];
	char str1[100];
	lint v_num=0, e_num=0, offset=0, g_v_size=0;
	lint v_id=0, v_to=0, e_weight=0;

	//open the file to read global information
	v_num=global_info[num_nodes];
	e_num=global_info[num_nodes+1];
	v_num_send = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	v_num_recv = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	vg_num_send = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	vg_num_recv = (lint*)_mm_malloc(sizeof(lint)*v_num, 64);
	for(i=0;i<v_num;i++){
		v_num_send[i] = 0;
		v_num_recv[i] = 0;
		vg_num_send[i] = 0;
		vg_num_recv[i] = 0;
	}

	//init mpi
	lint num_read =0;
	lint total = global_info[myrank];
	//init my own read file
	char extension[100];
	char extension1[100];
	snprintf(extension, 10, "%d", myrank);
	snprintf(extension1, 10, "%d", num_nodes);
	concate(f_buf, 7, dir, data_extend, file, "_phy_", extension , "_", extension1);
	if((reader=fopen(f_buf, "r"))==NULL){
		DPRINTF(1, "the file %s does not exists!\n", f_buf);
		ERROR_PRINT();
	}
	//read files
	while(!feof(reader)){
		fscanf(reader, "%s %lld %lld %lld", str, &v_id , &v_to, &e_weight);
		v_id -= 1;
		v_to -= 1;
		v_num_send[v_id]++;	
		v_num_recv[v_to]++;	
		num_read++;
		if(num_read == total)
			break;
	}
	//gather all the results
	MPI_Allreduce(v_num_send,vg_num_send,v_num,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);		
	MPI_Allreduce(v_num_recv,vg_num_recv,v_num,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);		
	/*FILE *writer = fopen("./mpi_vg_num", "w");
	for(i=0;i<v_num;i++)
		fprintf(writer, "%lld\n", vg_num_send[i]);
	fclose(writer);*/
	//compute distrib
	lint distrib_num = e_num/num_nodes;
	for(i=0;i<num_nodes;i++)
	{
		distrib[i]=0;
		send_size[i]=0;
		recv_size[i]=0;
	}
	for (i=0;i<v_num;i++)
	{
		send_size[j] += vg_num_send[i];
		recv_size[j] += vg_num_recv[i];
		if((send_size[j]>distrib_num && type==SEND) 
				|| (recv_size[j]>distrib_num && type == RECV))
		{
			distrib[j]=i;
			j++;
		}
	}
	distrib[num_nodes-1]=v_num-1;
	//write recip in parallel
	FILE *recip;	
	snprintf(extension, 10, "%d", myrank);
	concate(f_buf, 7, dir, par_extend, file, "_recip_", extension, "_", extension1);	
	if((recip=fopen(f_buf, "w"))==NULL){
		DPRINTF(1, "the file %s does not exists!\n", f_buf);
		ERROR_PRINT();
	}
	lint start = (myrank==0?0:distrib[myrank-1]);
	lint end = distrib[myrank];
//	DPRINTF(3, "e_num %d distrib_num %d start %d end %d\n", e_num, distrib_num, start, end);
	if(myrank==0){
		for(i=0;i<g_v_num;i++){
			DPRINTF(1, "%lld ", vg_num_send[i]);
		}
		DPRINTF(1, "\n");
	}
	printf("myrank %d start %lld end %lld\n", myrank, start, end);
	for(i=(start==0?0:(start+1)); i<=end; i++){
		fprintf(recip, "%lf ", 1.0/(double)vg_num_send[i]);
	}
	fprintf(recip, "\n");
	for(i=0;i<num_nodes;i++)
		fprintf(recip, "%lld ", distrib[i]);
	fprintf(recip, "\n\n");
	fclose(recip);
	_mm_free(v_num_send);	
	_mm_free(v_num_recv);	
	_mm_free(vg_num_send);	
	_mm_free(vg_num_recv);	
}

/**
 * @brief parition the graph to make sure that each subgraph has 
 * even number of outgoing edges.  
 *
 * compute the distribution of each node which has differnet number of 
 * vertices and euqal number of out going egdes
 *
 * @param file base file name
 * @param dir base directory
 * @param distrib how vertices are distribute among different files
 * @param send_size how many outgoing edges a subset of vertices has
 * @param recv_size how many ingoing edges a subset of vertices has
 * @param partition by #outgoing edges =>SEND #ingoing edges =>RECV
 */
	void 
chunk_files_par(char* file, char *dir, 
		lint *distrib, lint *send_size, 
		lint *recv_size, sint type, 
		lint *global_info, sint myrank)
{
	lint v_id, v_to, e_weight;
	lint i=0,j=0;
	sint bucket;
	lint offset =1;
	char f_buf[100];
	//write upper_bounds
	sprintf(f_buf, "%s%s%s_send_%d", dir, par_extend, file, num_nodes);
	lint **send_ub = (lint**)malloc(sizeof(lint*)*num_nodes);
	lint **g_send_ub = (lint**)malloc(sizeof(lint*)*num_nodes);
	for(i=0;i<num_nodes;i++){
		send_ub[i] = (lint*)malloc(sizeof(lint)*num_nodes);
		g_send_ub[i] = (lint*)malloc(sizeof(lint)*num_nodes);
		for(j=0;j<num_nodes;j++)
			send_ub[i][j]=0;
	}
	FILE *ub_writer;
	if((ub_writer = fopen(f_buf, "w"))==NULL){
		DPRINTF(1, "the file %s does not exist\n", f_buf);
		ERROR_PRINT();
	}
	//init mpi and know my rank
	int *argc = (int*)malloc(sizeof(int*));
	char ***argv = (char***)malloc(sizeof(char***));
	////////////////////////part one///////////////////////////////////////
	//read my own file and sort them into buckets
	FILE *reader;
	FILE *writer[num_nodes];	
	char str[100];
	char str1[100];
	char str2[100];
	snprintf(str, 10, "%d", myrank);
	snprintf(str2, 10, "%d", num_nodes);
	concate(f_buf, 7, dir, data_extend, file, "_phy_", str, "_", str2);
	DPRINTF(3, "I am %d, I am reading %s\n", myrank, f_buf);
	if((reader = fopen(f_buf, "r"))==NULL){
		DPRINTF(1, "the file %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}
	for(i=0;i<num_nodes; i++){
		snprintf(str1, 10, "%lld", i);
		concate(f_buf, 9, dir, data_extend, file, "_", str, "_", str1, "_",str2);	
		if((writer[i] = fopen(f_buf, "w"))==NULL){
			DPRINTF(1, "the file %s does not exist!\n", f_buf);
			ERROR_PRINT();
		}
	}
	//body computation
	lint total = global_info[myrank];
	lint num_read =0;
	while(!feof(reader)){
		fscanf(reader, "%s %lld %lld %lld", str, &v_id , &v_to, &e_weight);
		v_id -= 1;
		v_to -= offset;
		sint from_node = bin_search(distrib, v_id, 0,  num_nodes-1);
                sint to_node = bin_search(distrib, v_to, 0,  num_nodes-1);
                send_ub[from_node][to_node]+=1;
		num_read++;
		if(type == RECV)
			bucket = bin_search(distrib, v_to, 0, num_nodes-1);
		else if(type == SEND)
			bucket = bin_search(distrib, v_id, 0, num_nodes-1);
		fprintf(writer[bucket], "%s %lld %lld %lld\n", str, v_id+1 , v_to+offset, e_weight);
		if(num_read==total)
			break;
	}
	//close files
	fclose(reader);
	for(i=0;i<num_nodes; i++){
		fprintf(writer[i], "t 1 1 1\n");
		fclose(writer[i]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////part two///////////////////////////////////////
	//read buckets and merge them into different file
	FILE *buckreader[num_nodes];
	FILE *finalwriter;	
	snprintf(str, 10, "%d", myrank);
	snprintf(str1, 10, "%d", num_nodes);
	if(type==SEND)
		concate(f_buf, 7, dir, data_extend, file, 
				"_sd_", str , "_", str1);
	else if(type == RECV)
		concate(f_buf, 7, dir, data_extend, file, 
				"_rc_", str , "_", str1);

	DPRINTF(3, "I am %d, I am writing %s\n", myrank, f_buf);
	if((finalwriter = fopen(f_buf, "w"))==NULL){
		DPRINTF(1, "the file %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}
	for(i=0;i<num_nodes; i++){
		snprintf(str1, 10, "%lld", i);
		concate(f_buf, 9, dir, data_extend, file, "_", str1, "_", str, "_", str2);	
		if((buckreader[i] = fopen(f_buf, "r"))==NULL){
			DPRINTF(1, "the file %s does not exist!\n", f_buf);
			ERROR_PRINT();
		}
	}
	//body computation
	//should print headers for each file
	lint num_v = (myrank==0?distrib[0]+1:distrib[myrank]-distrib[myrank-1]);
	offset = (myrank==0?1:distrib[myrank-1]+2);
	lint global_v = global_info[num_nodes];
	lint num_e;
	if(type==SEND)
		num_e = send_size[myrank];
	else if(type == RECV)
		num_e = recv_size[myrank];
	fprintf(finalwriter, "p sp %lld %lld %lld %lld\n", num_v, num_e, offset, global_v);
	for(i=0;i<num_nodes;i++){
		while(!feof(buckreader[i])){
			fscanf(buckreader[i], "%s %lld %lld %lld", str, &v_id , &v_to, &e_weight);
			
			if(str[0]=='t')
				break;
			else
				fprintf(finalwriter, "%s %lld %lld %lld\n", str, v_id , v_to, e_weight);
		}
	}
	//gather all the results
	for(i=0;i<num_nodes;i++){
		MPI_Allreduce(send_ub[i],g_send_ub[i],num_nodes,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
		for(j=0;j<num_nodes;j++)
			fprintf(ub_writer, "%lld ", g_send_ub[i][j]);
		fprintf(ub_writer, "\n");
	}
	fclose(ub_writer);	
	//close files
	fclose(finalwriter);
	for(i=0;i<num_nodes; i++)
		fclose(buckreader[i]);
}

void remove_files(char* file, char* dir, sint myrank){
	char str[100];
	char str1[100];
	char str2[100];
	char f_buf[100];
	lint i;
	//delete all the files generated in part one
	snprintf(str, 10, "%d", myrank);
	snprintf(str2, 10, "%d", num_nodes);
	concate(f_buf, 7, dir, data_extend, file, "_phy_", str, "_", str2);
	DPRINTF(3, "deleting file %s\n", f_buf);
	if(remove(f_buf)!=0){
		DPRINTF(1, "error deleting file %s\n", f_buf);
		ERROR_PRINT();
	}
	for(i=0;i<num_nodes; i++){
		snprintf(str1, 10, "%lld", i);
		concate(f_buf, 9, dir, data_extend, file, "_", str, "_", str1, "_", str2);	
		DPRINTF(3, "deleting file %s\n", f_buf);
		if(remove(f_buf)!=0){
			DPRINTF(1, "error deleting file %s\n", f_buf);
			ERROR_PRINT();
		}
	}
}
