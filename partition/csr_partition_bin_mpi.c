#include "graph.h"
#include "csr_partition_bin_mpi.h"
#include "hashmap.h"

#undef OFFSET
#define OFFSET 0

void check_msg(lint *msg, int size, int num, int myrank, int ii){
	int i;
	for(i=0;i<size; i++){
		lint from = msg[i*2];
		if(from==num){
			printf("this is my rank %d and send to %d\n", myrank, ii);
			break;
		}
	}
}

void prin_msg(lint *msg, lint size, int type, int ii){
	FILE *f;
	char f_buf[100];
	sprintf(f_buf, "./debug/%d_%d", type, ii);
	if((f=fopen(f_buf, "w"))==NULL)
		printf("file does not exist\n");
	lint i;
	for(i=0;i<size;i++)
		fprintf(f, "%lld %lld\n", msg[i*2], msg[i*2+1]);
	fprintf(f, "\n");
}

void
remove_duplicate(char* file, char*dir, sint myrank){
	char f_buf[100];
	sprintf(f_buf, "%s%s_%d", dir, file, myrank);
	FILE *reader;
	if((reader=fopen(f_buf, "rb"))==NULL){
		DPRINTF(1, "THE FILE %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}
	FILE *writer;
	sprintf(f_buf, "%s%s%s_phy_%d_%d", dir, data_extend, file, myrank, num_nodes);
	if((writer=fopen(f_buf, "wb"))==NULL){
		DPRINTF(1, "THE FILE %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}

	lint edge_num;
	fread(&edge_num, sizeof(lint), 1, reader);
	lint *buf = (lint*)_mm_malloc(sizeof(lint)*2*edge_num, 64);
	sint *dup = (sint*)_mm_malloc(sizeof(sint)*edge_num, 64);
	fread(buf, sizeof(lint), 2*edge_num, reader);
	lint i,j;
	for(i=0;i<edge_num;i++){
		dup[i] = FALSE;
	}
	p_elem *hash_map = (p_elem*)malloc(sizeof(p_elem)*HASH_BASE);
	init_hash(HASH_BASE, hash_map);
	//major computaiton
	lint cur = buf[0]; //first edge
	lint start = 0;
	lint nondup_edge_num=0;
	for(i=0;i<edge_num;i++){
		lint from = buf[2*i];
		lint to = buf[2*i+1];
		if( (from == -1) || (to == -1)) { printf("Node: %d ::: i: %lld, graph500 generated graph has -1 in it", myrank, i); ERROR_PRINT(); }
		lint found = search_hash(from, to, hash_map);
		if(found == FALSE){
			dup[i]=FALSE;
			nondup_edge_num++;
			insert_hash(from, to, hash_map);
		}
		else{
			dup[i]=TRUE;
			//printf("Previous was at edge #: %lld / %lld on Node: %d\n", i, edge_num, myrank);
		}
	}
	free_hash(HASH_BASE, hash_map);
	//out put back
	lint off =1;
	if(num_nodes == 1){
		fwrite(&g_v_num, sizeof(lint), 1, writer);
	}
	fwrite(&nondup_edge_num, sizeof(lint), 1, writer);	
	if(num_nodes == 1){
		lint off = OFFSET;
		fwrite(&off, sizeof(lint), 1, writer);
		fwrite(&g_v_num, sizeof(lint), 1, writer);
	}
	lint *nondup_buf = (lint*)_mm_malloc(sizeof(lint)*2*nondup_edge_num, 64);
	j=0;
	for(i=0;i<edge_num;i++){
		if(dup[i]==FALSE){
			nondup_buf[2*j] = buf[2*i];
			nondup_buf[2*j+1] = buf[2*i+1];
			if((nondup_buf[2*j] == -1) || (nondup_buf[2*j+1] == -1)) ERROR_PRINT();
			j++;
		}
	}
	if(j != nondup_edge_num) ERROR_PRINT();
	fwrite(nondup_buf, sizeof(lint), nondup_edge_num*2, writer);
	fclose(reader);	
	fclose(writer);	
	_mm_free(buf);
	_mm_free(dup);
	_mm_free(nondup_buf);
}
void myMPI_Allreduce(lint *sendbuf, lint *recvbuf, int count, 
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
#if 0
	const long long int MAX_COUNT = (1<<25);
	int rounds = count / MAX_COUNT;
	if(rounds * MAX_COUNT != count) rounds++;
	int round;
	for(round = 0; round < rounds; round++)
	{
		int count_per_round = (round == (rounds-1))?(count - round*MAX_COUNT):(MAX_COUNT);
		MPI_Allreduce(sendbuf + round*count_per_round, recvbuf + round*count_per_round, count_per_round, datatype, op, comm);
	}
#else
	int rounds = 8;
	int count_per_round = count/rounds;
	if(count_per_round * rounds != count) ERROR_PRINT();
	int round;
	for(round = 0; round < rounds; round++)
	{	
		MPI_Allreduce(sendbuf + round*count_per_round, recvbuf + round*count_per_round, count_per_round, datatype, op, comm);
	}
#endif
}
void comp_distrib_csr_par_bin(char* file, char* dir, 
		lint *distrib, lint *send_size,
		lint *recv_size, sint type, sint myrank){
	lint i=0,j=0;
	lint *v_num_send;
	lint *v_num_recv;
	lint *vg_num_send;
	lint *vg_num_recv;
	FILE *reader;	
	char f_buf[100];
	lint v_num=0, e_num=0, offset=0;
	lint v_id=0, v_to=0, e_weight=0;

	//open the file to read global information

	v_num_send = (lint*)_mm_malloc(sizeof(lint)*g_v_num, 64);
	v_num_recv = (lint*)_mm_malloc(sizeof(lint)*g_v_num, 64);
	vg_num_send = (lint*)_mm_malloc(sizeof(lint)*g_v_num, 64);
	vg_num_recv = (lint*)_mm_malloc(sizeof(lint)*g_v_num, 64);
	for(i=0;i<g_v_num;i++){
		v_num_send[i] = 0;
		v_num_recv[i] = 0;
		vg_num_send[i] = 0;
		vg_num_recv[i] = 0;
	}

	//init mpi
	sprintf(f_buf, "%s%s%s_phy_%d_%d", dir, data_extend, file, myrank, num_nodes);
	if((reader=fopen(f_buf, "rb"))==NULL){
		DPRINTF(1, "the file %s does not exists!\n", f_buf);
		ERROR_PRINT();
	}
	//read files
	fread(&e_num, sizeof(lint), 1, reader);
	lint *edge_arr = (lint*)_mm_malloc(sizeof(lint)*e_num*2, 64); 
	fread(edge_arr, sizeof(lint), e_num*2, reader);
	for(i=0;i<e_num;i++){
		lint v_from = edge_arr[2*i]-OFFSET;
		lint v_to = edge_arr[2*i+1]-OFFSET;
		v_num_send[v_from]++;
		v_num_recv[v_to]++;
	}
	//gather all the results
	if(v_num_send == NULL) ERROR_PRINT();
	if(vg_num_send == NULL) ERROR_PRINT();
	myMPI_Allreduce(v_num_send,vg_num_send,g_v_num,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
	myMPI_Allreduce(v_num_recv,vg_num_recv,g_v_num,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);		
	/*FILE *writer = fopen("./bin_vg_num", "w");
	for(i=0;i<g_v_num;i++)
		fprintf(writer, "%lld\n", vg_num_send[i]);
	fclose(writer);*/
	//compute distrib
	lint total_e_num=0;
	for(i=0;i<g_v_num;i++){
		total_e_num+=vg_num_send[i];
	}
	lint distrib_num = total_e_num/num_nodes;
	//printf("total e_num is %lld")
	for(i=0;i<num_nodes;i++)
	{
		distrib[i]=0;
		send_size[i]=0;
		recv_size[i]=0;
	}
	for (i=0;i<g_v_num;i++)
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
	distrib[num_nodes-1]=g_v_num-1;
	//if(myrank==0)
	//{
	//	int i;
	//	for(i = 0; i < num_nodes; i++) DPRINTF(1, "distrib[%d] = %lld\n", i, distrib[i]);
	//	//DPRINTF(1, "the last distrib is: %lld\n", distrib[num_nodes-1]);
	//}
	////write recip in parallel
	FILE *recip;	
	sprintf(f_buf, "%s%s%s_recip_%d_%d", dir, par_extend, file, myrank, num_nodes);
	if((recip=fopen(f_buf, "wb"))==NULL){
		DPRINTF(1, "the file %s does not exists!\n", f_buf);
		ERROR_PRINT();
	}

	lint start = (myrank==0?0:distrib[myrank-1]);
	lint end = distrib[myrank];
	lint size=end-start+(start==0?1:0);
	//printf("Node: %d asking for %lld bytes\n", myrank, (sizeof(double)*size));
	double *recip_buf = (double*)_mm_malloc(sizeof(double)*size, 64);
	for(i=0,j=(start==0?0:(start+1)); j<=end; i++,j++){
		recip_buf[i]=1.0/vg_num_send[j];
	}
	//printf("myrank %d start %lld end %lld size %lld last elem %lld\n", myrank, start, end, size, g_num_send[end]);
	//DPRINTF(1, "LAST ELEM IS %f %lld\n", recip_buf[end], end);
	//if(myrank==0)
	//for(i=0;i<g_v_num;i++)
	//	DPRINTF(1, "%lld|%f ", i, recip_buf[i]);
	//DPRINTF(1, "\n");
	fwrite(distrib, sizeof(lint), num_nodes, recip);
	fwrite(recip_buf, sizeof(double), size, recip);
	fclose(recip);
#ifdef DEBUG_BIN
	printf("offset is %d\n", OFFSET);
	sprintf(f_buf, "./bin2t/parex/%s_%d", "recip", myrank);
	FILE *brecip = fopen(f_buf, "w");
	for(i=0; i<size; i++)
		fprintf(brecip, "%lf ", recip_buf[i]);
	fprintf(brecip, "\n");
	for(i=0;i<num_nodes;i++)
		fprintf(brecip, "%lld ", distrib[i]);
	fprintf(brecip, "\n");
	fclose(brecip);	
#endif	
	

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
	if((ub_writer = fopen(f_buf, "wb"))==NULL){
		DPRINTF(1, "the file %s does not exist\n", f_buf);
		ERROR_PRINT();
	}
	for(i=0;i<e_num;i++){
		lint v_from = edge_arr[2*i]-OFFSET;
		lint v_to = edge_arr[2*i+1]-OFFSET;
		sint from_node = bin_search(distrib, v_from, 0,  num_nodes-1);
		sint to_node = bin_search(distrib, v_to, 0,  num_nodes-1);
		send_ub[from_node][to_node]+=1;
	}
	//gather all the results
	for(i=0;i<num_nodes;i++){
		MPI_Allreduce(send_ub[i],g_send_ub[i],num_nodes,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
		fwrite(g_send_ub[i], sizeof(lint), num_nodes, ub_writer);
	}
#ifdef DEBUG_BIN
	sprintf(f_buf, "./bin2t/parex/ub_%d", num_nodes);
	FILE *bub = fopen(f_buf, "w");
	for(i=0;i<num_nodes;i++){
		for(j=0; j<num_nodes; j++)
			fprintf(bub, "%lld ", g_send_ub[i][j]);
		fprintf(bub, "\n");
	}
	fclose(bub);
#endif
	fclose(ub_writer);

	//_mm_free(v_num_send);	
	//_mm_free(v_num_recv);	
	//_mm_free(vg_num_send);	
	//_mm_free(vg_num_recv);	
	//_mm_free(edge_arr);	
}

void read_send_ub_bin(char *file, char *dir, lint **send_ub){
	lint i,j;
	FILE *f;
	char buf[100];
	sprintf(buf, "%s%s%s_send_%s", dir, par_extend, file, nodes_extension);
	if((f=fopen(buf, "rb"))==NULL){
		printf("the file %s does not exists!\n", buf);
		ERROR_PRINT();
	}
	for(i=0;i<num_nodes;i++)
		fread(send_ub[i], sizeof(lint), num_nodes, f);	
	fclose(f);
}

void read_recip_bin(char* file, char *dir, lint *distrib, 
		t_csr *graph, sint myrank){
	lint i,j,k;
	FILE *reader;
	//FILE *writer;
	char f_buf[100];
	double in;
	sprintf(f_buf, "%s%s%s_recip_%d_%d", dir, par_extend, file, myrank, num_nodes);
	if((reader = fopen(f_buf, "rb"))==NULL){
		printf("the file does not exists!\n");
		ERROR_PRINT();
	}
	double *recip_read = (double*)_mm_malloc(sizeof(double)*graph->v_size, 64);
	fread(distrib, sizeof(lint), num_nodes, reader);
	fread(recip_read, sizeof(double), graph->v_size, reader);
	for(i=0;i<graph->v_size;i++){
		graph->vet_info[i].recip = recip_read[i];
	}
	_mm_free(recip_read);
	fclose(reader);
}


//	void 
//chunk_files_bin_par_file(char* file, char *dir, 
//		lint *distrib, sint type, 
//		sint myrank){
//	lint v_id, v_to, e_weight;
//	lint i=0,j=0;
//	sint bucket;
//	lint offset =1;
//	////////////////////////part one///////////////////////////////////////
//	//read my own file and sort them into buckets
//	if(myrank == 0) printf("enter here\n");
//	FILE *reader;
//	char f_buf[100];
//	sprintf(f_buf, "%s%s%s_phy_%d_%d", dir, data_extend, file, myrank, num_nodes);
//	if((reader = fopen(f_buf, "rb"))==NULL){
//		DPRINTF(1, "the file %s does not exist!\n", f_buf);
//		ERROR_PRINT();
//	}
//	//body computation
//	int v_num, e_num;
//	fread(&e_num, sizeof(lint), 1, reader);
//	lint *edge_arr = (lint*)_mm_malloc(sizeof(lint)*e_num*2, 64); 
//	lint **buffer = (lint**)malloc(sizeof(lint*)*num_nodes);	
//	lint *idx = (lint*)malloc(sizeof(lint)*num_nodes);
//	for(i=0;i<num_nodes;i++){
//		idx[i]=0;
//		buffer[i] = (lint*)malloc(sizeof(lint)*e_num*2);
//	}
//	fread(edge_arr, sizeof(lint), e_num*2, reader);
//	for(i=0;i<e_num;i++){
//		lint v_from = edge_arr[2*i]-OFFSET;
//		lint v_to = edge_arr[2*i+1]-OFFSET;
//		//if( (v_from) < 0 || (v_to < 0)) ERROR_PRINT();
//		if(type == RECV){
//			bucket = bin_search(distrib, v_to, 0, num_nodes-1);
//			int lower_bound = (bucket == 0)?0:distrib[bucket-1];
//			int upper_bound = distrib[bucket];
//			//if( (v_to < lower_bound) ||  (v_to > upper_bound)) ERROR_PRINT();
//		}
//		else if(type == SEND){
//			bucket = bin_search(distrib, v_from, 0, num_nodes-1);
//			lint lower_bound = (bucket == 0)?0:distrib[bucket-1];
//			lint upper_bound = distrib[bucket];
//		}
//		buffer[bucket][idx[bucket]*2]=v_from;
//		buffer[bucket][idx[bucket]*2+1]=v_to;
//		idx[bucket]++;
//	}
//	FILE *writer[num_nodes];
//	for(i=0;i<num_nodes;i++){
//		sprintf(f_buf, "%s%s%s_%d_%lld_%d", dir, data_extend, file, myrank, i, num_nodes);
//		if((writer[i] = fopen(f_buf, "wb"))==NULL){
//			DPRINTF(1, "the file %s does not exist!\n", f_buf);
//			ERROR_PRINT();
//		}
//		fwrite(&idx[i], sizeof(lint), 1, writer[i]);
//		fwrite(buffer[i], sizeof(lint), idx[i]*2, writer[i]);
//		fclose(writer[i]);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	FILE *reader1[num_nodes];
//	FILE *finalwriter;
//	if(type==SEND){
//		sprintf(f_buf, "%s%s%s_sd_%d_%d", dir, data_extend, file, myrank, num_nodes);
//	}
//	else if(type ==RECV){
//		sprintf(f_buf, "%s%s%s_rc_%d_%d", dir, data_extend, file, myrank, num_nodes);
//	}
//	if((finalwriter = fopen(f_buf, "wb"))==NULL){
//		printf("the file %s does not exists\n", f_buf);
//		ERROR_PRINT();
//	}
//	lint total_e_num=0;
//	lint *e_n = (lint*)malloc(sizeof(lint)*num_nodes);
//	for(i=0;i<num_nodes;i++){
//		sprintf(f_buf, "%s%s%s_%lld_%d_%d", dir, data_extend, file, i, myrank,num_nodes);
//		if((reader1[i] = fopen(f_buf, "rb"))==NULL){
//			DPRINTF(1, "the file %s does not exist!\n", f_buf);
//			ERROR_PRINT();
//		}
//		int start =0;
//		fread(&e_n[i], sizeof(lint), 1, reader1[i]);
//		total_e_num+=e_n[i];
//	}
//	lint *e = (lint*)malloc(sizeof(lint)*total_e_num*2);
//	lint start =0;
//	for(i=0;i<num_nodes;i++){
//		fread(e+(start*2), sizeof(lint), e_n[i]*2, reader1[i]);
//		start+=e_n[i];
//		fclose(reader1[i]);
//	}
//	lint num_v = (myrank==0?distrib[0]+1:distrib[myrank]-distrib[myrank-1]);
//	offset = (myrank==0?0:distrib[myrank-1]+1)+OFFSET;
//	fwrite(&num_v, sizeof(lint), 1, finalwriter);
//	fwrite(&total_e_num, sizeof(lint), 1, finalwriter);
//	fwrite(&offset, sizeof(lint), 1, finalwriter);
//	fwrite(&g_v_num, sizeof(lint), 1, finalwriter);
//	printf("%d writing %lld %lld %lld %lld\n", myrank, num_v, total_e_num, offset, g_v_num);
//	fflush(stdout);
//	fwrite(e, sizeof(lint), total_e_num*2, finalwriter);
//	fclose(finalwriter);
//}
	void 
chunk_files_bin_par(char* file, char *dir, 
		lint *distrib, sint type, 
		sint myrank)
{
	lint v_id, v_to, e_weight;
	lint i=0,j=0;
	sint bucket;
	lint offset =1;
	////////////////////////part one///////////////////////////////////////
	//read my own file and sort them into buckets
	FILE *reader;
	FILE *writer[num_nodes];	
	char f_buf[100];
	sprintf(f_buf, "%s%s%s_phy_%d_%d", dir, data_extend, file, myrank, num_nodes);
	if((reader = fopen(f_buf, "rb"))==NULL){
		DPRINTF(1, "the file %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}
	//body computation
	int v_num, e_num;
	fread(&e_num, sizeof(lint), 1, reader);
	lint *edge_arr = (lint*)_mm_malloc(sizeof(lint)*e_num*2, 64); 
	sint *node_mark = (sint*)_mm_malloc(sizeof(sint)*e_num, 64); 
	lint *msg = (lint*)_mm_malloc(sizeof(lint)*e_num*2, 64); 
	fread(edge_arr, sizeof(lint), e_num*2, reader);
	//DPRINTF(1, "num edges %lld and the first edge is %lld %lld\n", e_num,
	//edge_arr[2*i]-1, edge_arr[2*i+1]-1);
	//DPRINTF(1, "%d %d %d %d\n", distrib[0], distrib[1],distrib[2],distrib[3]);
	//printf("OFFSET: %d\n", OFFSET);
	for(i=0;i<e_num;i++){
		lint v_from = edge_arr[2*i]-OFFSET;
		lint v_to = edge_arr[2*i+1]-OFFSET;
		//if( (v_from) < 0 || (v_to < 0)) ERROR_PRINT();
		if(type == RECV){
			bucket = bin_search(distrib, v_to, 0, num_nodes-1);
			int lower_bound = (bucket == 0)?0:distrib[bucket-1];
			int upper_bound = distrib[bucket];
			//if( (v_to < lower_bound) ||  (v_to > upper_bound)) ERROR_PRINT();
		}
		else if(type == SEND){
			bucket = bin_search(distrib, v_from, 0, num_nodes-1);
			lint lower_bound = (bucket == 0)?0:distrib[bucket-1];
			lint upper_bound = distrib[bucket];
			/*if( (v_from < lower_bound) ||  (v_from > upper_bound)) 
			  {
			  printf("Node: %d ::: i: %lld v_from: %lld, bucket: %d, bucket_lower_bound: %lld, bucket_upper_bound: %lld\n", myrank, i, v_from, bucket, lower_bound, upper_bound);
			  ERROR_PRINT();
			  }*/
		}
		node_mark[i]=bucket;	
		//if(v_from == 260662193 && myrank ==12)
		//	DPRINTF(1, "%d %d\n", bucket, myrank);
		//DPRINTF(1, "%d %d %d\n", v_from, v_to, bucket);
	}
	////////////////////////part two///////////////////////////////////////
	//check msg size
	MPI_Request *request = (MPI_Request*)malloc(sizeof(MPI_Request)*num_nodes);
	lint *msg_size = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
	lint **msg_recv_size = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	for(i=0;i<num_nodes;i++)
		msg_recv_size[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
	lint total_recv_size=0;
	bin_check_msg_size(msg_size, node_mark, e_num);
	for(i=0; i<num_nodes; i++){
		if(i==myrank)
			continue;

		//MPI_Isend(msg_size, num_nodes, MPI_LONG, i, 4*num_nodes + myrank, MPI_COMM_WORLD, &request);	
		MPI_Isend(msg_size, num_nodes, MPI_LONG, i, 4, MPI_COMM_WORLD, &request[i]);	
	}
	for(i=0; i<num_nodes; i++){
		if(i==myrank)
			continue;
		//MPI_Recv(msg_recv_size[i], num_nodes, MPI_LONG, i, 4*num_nodes + i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		MPI_Recv(msg_recv_size[i], num_nodes, MPI_LONG, i, 4, MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
		total_recv_size+=msg_recv_size[i][myrank];
	}
	total_recv_size+=msg_size[myrank];
	lint *edge_recv_arr = (lint*)_mm_malloc(sizeof(lint)*total_recv_size*2, 64); 
	//DPRINTF(1, "total num of edges in %d is %lld\n", myrank, total_recv_size);
	int num_iter = 64;	
	int ii=0;
	lint **send_size_div = (lint**)malloc(sizeof(lint*)*num_nodes);
	lint **recv_size_div = (lint**)malloc(sizeof(lint*)*num_nodes);
	for(i=0;i<num_nodes;i++){
		send_size_div[i]=(lint*)malloc(sizeof(lint)*num_iter);
		recv_size_div[i]=(lint*)malloc(sizeof(lint)*num_iter);
		lint avg_send = msg_size[i]/num_iter+1;
		lint avg_recv = msg_recv_size[i][myrank]/num_iter+1;
		lint remain_send= msg_size[i];
		lint remain_recv= msg_recv_size[i][myrank];
		for(j=0;j<num_iter;j++){
			send_size_div[i][j]=remain_send>avg_send?avg_send:remain_send;	
			recv_size_div[i][j]=remain_recv>avg_recv?avg_recv:remain_recv;	
			remain_send-=avg_send;
			remain_recv-=avg_recv;
		}
	}
	//printf("the size is %lld\n", send_size_div[0][0]);
	lint *send_start = (malloc(sizeof(lint)*num_nodes));
	for(i=0;i<num_nodes;i++)
		send_start[i]=0;
	lint start =0;
	for(ii=0;ii<num_iter;ii++){
		for(i=0; i<num_nodes; i++){
			if(i==myrank)
				continue;
			send_start[i] = bin_to_msg(edge_arr, node_mark, msg, e_num, i, send_start[i], send_size_div[i][ii]);
			//check_msg(msg, msg_recv_size[i][myrank], 260662193, myrank, i);
			//MPI_Isend(msg, size*2, MPI_LONG, i, 5*num_nodes + myrank, MPI_COMM_WORLD, &request);	
			MPI_Send(msg, send_size_div[i][ii]*2, MPI_LONG_LONG, i, 5, MPI_COMM_WORLD);
		}	
		for(i=0;i<num_nodes; i++){
			if(i==myrank){
				send_start[i] = bin_to_msg(edge_arr, node_mark, msg, e_num, i, send_start[i], send_size_div[i][ii]);
				start = bin_from_msg(edge_recv_arr, msg, send_size_div[i][ii], start);
				continue;
			}

			//MPI_Recv(msg, msg_recv_size[i][myrank]*2, MPI_LONG, i, 5*num_nodes + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			MPI_Recv(msg, recv_size_div[i][ii]*2, MPI_LONG_LONG, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
			start = bin_from_msg(edge_recv_arr, msg, recv_size_div[i][ii], start);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	FILE *finalwriter;	

	if(type==SEND)
		sprintf(f_buf, "%s%s%s_sd_%d_%d", dir, data_extend, file, myrank, num_nodes);
	else if(type == RECV)
		sprintf(f_buf, "%s%s%s_rc_%d_%d", dir, data_extend, file, myrank, num_nodes);
	if((finalwriter = fopen(f_buf, "wb"))==NULL){
		DPRINTF(1, "the file %s does not exist!\n", f_buf);
		ERROR_PRINT();
	}
	lint num_v = (myrank==0?distrib[0]+1:distrib[myrank]-distrib[myrank-1]);
	offset = (myrank==0?0:distrib[myrank-1]+1)+OFFSET;
	fwrite(&num_v, sizeof(lint), 1, finalwriter);
	fwrite(&total_recv_size, sizeof(lint), 1, finalwriter);
	fwrite(&offset, sizeof(lint), 1, finalwriter);
	fwrite(&g_v_num, sizeof(lint), 1, finalwriter);
	//for(i=0;i<10;i++)
	//	printf("%lld|%lld ", edge_recv_arr[2*i], edge_recv_arr[2*i+1]);
	//printf("\n");	
	fwrite(edge_recv_arr, sizeof(lint), total_recv_size*2, finalwriter);
	//DPRINTF(1, "finished writing files %lld %lld %lld %lld\n", num_v, total_recv_size, offset, g_v_num);
	//if(myrank==1){
	//	DPRINTF(1, "%lld %lld %lld %lld\n", num_v, total_recv_size, offset, g_v_num);
	//	for(i=0;i<total_recv_size; i++)
	//		DPRINTF(1, "a %lld %lld\n", edge_recv_arr[2*i], edge_recv_arr[2*i+1]);
	//}
#ifdef DEBUG_BIN
	lint *key = _mm_malloc(sizeof(lint)*total_recv_size, 64);
	lint *val = _mm_malloc(sizeof(lint)*total_recv_size, 64);
	for(i=0;i<total_recv_size;i++){
		key[i] = edge_recv_arr[2*i];
		val[i] = edge_recv_arr[2*i+1];
	}
	printf("start sorting, the size is %lld\n", total_recv_size);
	q_sort_two(key, val, 0, total_recv_size-1);
	//perform local sort
	start=0;
	lint end=0;
	while(1){
		while(1){
			if(key[end]!=key[start] || end==total_recv_size)
				break;
			end++;
		}
		q_sort_two(val, key, start, end-1);
		if(end == total_recv_size)
			break;
		start = end;
	}
	sprintf(f_buf, "./bin2t/datex/%s_%d", type==SEND?"sd":"rc", myrank);
	FILE *datf = fopen(f_buf, "w");
	fprintf(datf, "p sp %lld %lld %lld %lld\n", num_v, total_recv_size, offset, g_v_num);
	for(i=0; i<total_recv_size; i++)
		fprintf(datf, "%lld %lld\n", key[i]+1, val[i]+1);
	fclose(datf);
#endif
	printf("%d writing %lld %lld %lld %lld\n", myrank, num_v, total_recv_size, offset, g_v_num);
	fclose(finalwriter);
}

void bin_remove_files(char* file, char* dir, sint myrank){
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
}

int
bin_to_msg(lint *edge_arr, sint *node_mark, 
		lint *msg, lint len, 
		sint node_id, lint start, lint size){
	lint i, j=0;
	for(i=start; i<len; i++){
		if(j==size)
			return i;
		if(node_mark[i]==node_id){
			msg[j*2]=edge_arr[i*2];
			msg[j*2+1]=edge_arr[i*2+1];
			j++;
		}
	}
	return i;
}

int
bin_from_msg(lint *edge_recv_arr, lint *msg, 
		lint len, lint start){
	lint i=0,j=start;
	for(i=0;i<len;i++,j++){
		edge_recv_arr[j*2]=msg[i*2];	
		edge_recv_arr[j*2+1]=msg[i*2+1];	
	}
	return j;	
}

void 
bin_check_msg_size(lint *msg_size, sint *node_mark, 
		lint len){
	lint i;
	for(i=0;i<num_nodes;i++)
		msg_size[i]=0;
	for(i=0;i<len;i++){
		int id = node_mark[i];
		msg_size[id]++;
	}	
}
