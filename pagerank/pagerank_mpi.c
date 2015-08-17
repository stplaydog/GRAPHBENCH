#include "graph.h"
#include "pagerank_mpi.h"
#include "utils.h"
//#include "taskQ.h"
#include "timer.h"
#include "pthread.h"
#include <immintrin.h>

//sint ***dirty;
t_time_list *t_list;
double *diff_tmp;
lint *idx_tmp;

void 
print_result(p_time_list pl,  int num_iter, 
		int num_edges, int recv_size);
void
reduce_msg(t_msg *msg, sint myrank, 
		lint **send_ub);
void to_msg_no_buf(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank);

	void 
print_result(p_time_list pl,  int num_iter, 
		int num_edges, int recv_size)
{
	printf("==========LOCAL COMPUATION============\n");
	printf("cycles: %lld  cycles per iteration: %lld cycles per iteration per edge: %lld\n", pl->list[0].cycle, pl->list[0].cycle/num_iter, pl->list[0].cycle/num_iter/num_edges);
	printf("time: %lf  time per iteration: %lf time per iteration per edge: %lf\n", pl->list[0].time, pl->list[0].time/num_iter, pl->list[0].time/num_iter/num_edges);
	printf("bandwidth achieved in Gb: %lf\n", 83.0e-9/( pl->list[0].time/num_iter/num_edges));
	printf("==========TO MESSAGE============\n");
	printf("cycles: %lld  cycles per iteration: %lld cycles per iteration per edge: %lld\n", pl->list[1].cycle, pl->list[1].cycle/num_iter, pl->list[1].cycle/num_iter/num_edges);
	printf("time: %lf  time per iteration: %lf time per iteration per edge: %lf\n", pl->list[1].time, pl->list[1].time/num_iter, pl->list[1].time/num_iter/num_edges);
	printf("bandwidth achieved in Gb: %lf\n", 83.0e-9/( pl->list[1].time/num_iter/num_edges));
	printf("==========MESSAGE PASSING============\n");
	printf("cycles: %lld  cycles per iteration: %lld cycles per iteration per edge: %lld\n", pl->list[2].cycle, pl->list[2].cycle/num_iter, pl->list[2].cycle/num_iter/num_edges);
	printf("time: %lf  time per iteration: %lf time per iteration per edge: %lf\n", pl->list[2].time, pl->list[2].time/num_iter, pl->list[2].time/num_iter/(num_edges+recv_size));
	printf("mpi bandwidth achieved in Gb: %lf\n", 16e-9/( pl->list[2].time/num_iter/(num_edges+recv_size)));
	printf("==========FROM MESSAGE============\n");
	printf("cycles: %lld  cycles per iteration: %lld cycles per iteration per edge: %lld\n", pl->list[3].cycle, pl->list[3].cycle/num_iter, pl->list[3].cycle/num_iter/recv_size);
	printf("time: %lf  time per iteration: %lf time per iteration per edge: %lf\n", pl->list[3].time, pl->list[3].time/num_iter, pl->list[3].time/num_iter/recv_size);
	printf("bandwidth achieved in Gb: %lf\n", 83.0e-9/( pl->list[3].time/num_iter/recv_size));
	printf("==========CHECK TERMINATION============\n");
	printf("total time: %lf  time per iteration: %lf\n", pl->list[4].time, pl->list[4].time/num_iter);
	printf("cycles for pre to_msg %lld, cycles for post to_msg %lld\n", pl->list[5].cycle, pl->list[6].cycle);
}

void
init_t_buff(t_buff *buf, int size){
	sint i, j;
	for(i=0;i<num_threads;i++){
		buf->int_buffer[i] = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
		buf->double_buffer[i] = (double**)_mm_malloc(sizeof(double*)*num_nodes, 64);
		buf->buff_idx[i] = (sint*)_mm_malloc(sizeof(sint)*num_nodes, 64);
		for(j=0;j<num_nodes;j++){
			buf->buff_idx[i][j] = 0;
			buf->int_buffer[i][j] = (lint*)_mm_malloc(sizeof(lint)*size, 64);
			buf->double_buffer[i][j] = (double*)_mm_malloc(sizeof(double)*size, 64);
		}
	}
}

void 
init_weight(t_csr *gr){
	int i;
	for(i=0;i<gr->v_size;i++)
		gr->vet_info[i].weight = 1.0;
}


void init_msg(t_msg *msg, lint **send_ub, sint myrank){
	lint i, j;
	msg->double_send = (double**)malloc(sizeof(double*)*num_nodes);
	msg->double_recv = (double**)malloc(sizeof(double*)*num_nodes);
	msg->int_send = (lint**)malloc(sizeof(lint*)*num_nodes);
	msg->int_recv = (lint**)malloc(sizeof(lint*)*num_nodes);
	msg->idx_send = (lint*)malloc(sizeof(lint)*num_nodes);
	msg->idx_recv = (lint*)malloc(sizeof(lint)*num_nodes);
	for(i=0;i<num_nodes;i++){
		msg->double_send[i] = (double*)malloc(sizeof(double)*send_ub[myrank][i]);
		msg->double_recv[i] = (double*)malloc(sizeof(double)*send_ub[i][myrank]);
		msg->int_send[i] = (lint*)malloc(sizeof(lint)*send_ub[myrank][i]);
		msg->int_recv[i] = (lint*)malloc(sizeof(lint)*send_ub[i][myrank]);
		msg->idx_send[i]=0;
		msg->idx_recv[i]=0;
		for(j=0;j<send_ub[myrank][i];j++){
			msg->double_send[i][j]=0;
			msg->int_send[i][j]=0;
		}
		for(j=0;j<send_ub[i][myrank];j++){
			msg->double_recv[i][j]=0;
			msg->int_recv[i][j]=0;
		}
	}	
}

void to_msg_no_buf(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank){
	lint i,j,k;
	lint **pfix_sum = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);
	lint **pfix_idx = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);
	for(i=0; i<num_threads; i++){
		pfix_sum[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		pfix_idx[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		for(j=0;j<num_nodes;j++){
			pfix_sum[i][j]=0;
			pfix_idx[i][j]=0;
		}
	}
	//#pragma omp parallel shared(gs, msg, diff, pfix_sum, pfix_idx) private(i,j,k)
	//	{
	//sint tid = omp_get_thread_num();
	//printf("i am :%d\n", tid);
	tic_sr(t_list, 5);
#ifdef USE_OMP
	omp_set_num_threads(num_threads);
#pragma omp parallel for shared (gs, pfix_sum, myrank) 
#endif
	for(i=0; i<gs->v_size; i++){
		//printf("%d\n", i);
		sint tid= omp_get_thread_num();
		for(j=(i==0?0:gs->vet_idx[i-1]); j<gs->vet_idx[i]; j++){
			//there is no need to prefetch
			lint to = gs->edge_info[j].node_id;
			if(to != myrank)
				pfix_sum[tid][to]++;
		}
	}
#pragma omp barrier
	toc_sr(t_list, 5);
	//printf("cycles taken %lld\n", t_list->list[5].cycle);
	for(i=0; i<num_nodes; i++){
		lint sum = 0;
		pfix_idx[0][i]=0;
		for(j=1; j<num_threads; j++){
			sum += pfix_sum[j-1][i];
			pfix_idx[j][i]=sum;
			//if(myrank==0)
			//printf("%lld ", pfix_idx[i][j]);
		}
		//if(myrank==0)
		//printf("\n");
	}
	tic_sr(t_list, 6);
#ifdef USE_OMP
	omp_set_num_threads(num_threads);
#pragma omp parallel for shared (gs, pfix_idx) private(i, j)
#endif
	for(i=0; i<gs->v_size; i++){
		sint tid=omp_get_thread_num();
		lint source = i;
		double target_val = diff[source];
		lint from = myrank;
		for(j=(i==0?0:gs->vet_idx[i-1]); j<gs->vet_idx[i]; j++){
			lint target = gs->edge_idx[j];
			lint to = gs->edge_info[j].node_id;
			
			if(from!=to){
				lint pos = pfix_idx[tid][to];
				msg->int_send[to][pos]=target;
				msg->double_send[to][pos]=target_val;
				pfix_idx[tid][to]++;
			}	
		}
	}
	toc_sr(t_list, 6);

	//	}
	//#endif
}

void assign_diff(double *diff_tmp, lint *diff_assign,
			double *diff_pr, t_csr *gs, sint myrank)
{
	lint i, j;
#pragma omp parallel for num_threads(num_threads) shared(diff_tmp, diff_pr, gs) private(j)
	for(i=0; i<gs->v_size; i++){
		lint source = i;
		double target_val = diff_pr[source];
		for(j=(i==0?0:gs->vet_idx[i-1]); j<gs->vet_idx[i]; j++){
			lint assign_pos = diff_assign[j];
			diff_tmp[assign_pos] = target_val;
		}
	}
	//char f_buf[100];
	//sprintf(f_buf, "data/rand/log/2_20_%d", myrank);
	//FILE *log = fopen(f_buf, "w"); 
	//if(log==NULL)
	//	ERROR_PRINT();
	//for(i=0;i<gs->e_size; i++)
	//	fprintf(log, "%lf ", diff_tmp[i]);
	//fprintf(log, "\n");
}

// pidx_tmp, idx_tmp, node_idx_tmp are all partitioned once, we only need the index of diff_tmp
void to_msg_reduce(double *diff_tmp, t_msg *msg,
			lint *pidx_tmp, lint *idx_tmp, 
			lint *bin_idx, lint *node_idx_tmp){
	lint i,j;
#pragma omp parallel for num_threads(num_threads) private(j)
	for(i=0;i<HASH_BINS; i++){
		for(j=bin_idx[i];j<bin_idx[i+1];j++){
			lint node = node_idx_tmp[j];
			lint idx = idx_tmp[j];
			lint position = pidx_tmp[j];
			double value = diff_tmp[j];
			msg->int_send[node][position] = idx;
			msg->double_send[node][position] += value;
		}
	}
}

void set_msg(t_msg *msg, lint **send_ub, sint myrank){
	lint i, j;
#pragma omp parallel for num_threads(num_threads) private(i, j)
	for(i=0;i<num_nodes; i++){
		for(j=0;j<send_ub[myrank][i]; j++){
			msg->double_send[i][j]=0;
		}
	}
}

void to_msg_reduce_serial(double *diff_tmp, t_msg *msg,
			lint *pidx_tmp, lint *idx_tmp,
			lint *node_idx_tmp, lint size){
	lint i;
	for(i=0; i<size; i++){
		lint node = node_idx_tmp[i];
		lint idx = idx_tmp[i];
		lint position = pidx_tmp[i];
		double value = diff_tmp[i];
		msg->int_send[node][position] = idx;
		msg->double_send[node][position] += value;
	}

}

void reverse_diff_assign(lint *diff_assign, lint e_size){
	lint i;
	lint *tmp = (lint*)_mm_malloc(sizeof(lint)*e_size, 64);
	for(i=0;i<e_size; i++){
		lint pos = diff_assign[i];
		tmp[pos] = i;
	}
	for(i=0; i<e_size; i++)
		diff_assign[i] = tmp[i];
	_mm_free(tmp);
}

void partition_msg(lint *diff_assign, lint *pidx_tmp, 
			lint *idx_tmp, lint *bin_idx, 
			lint * node_idx_tmp, lint size){
	double *tmp_diff_assign = (double*)_mm_malloc(sizeof(double)*size, 64);	
	lint *tmp_pidx = (lint*)_mm_malloc(sizeof(lint)*size, 64);	
	lint *tmp_idx = (lint*)_mm_malloc(sizeof(lint)*size, 64);	
	lint *tmp_node_idx = (lint*)_mm_malloc(sizeof(lint)*size, 64);	

	lint N = size;
	sint *histogram = (sint*) _mm_malloc(sizeof(int)*NUMENTRIESPERTABLE*num_threads, 64);
	sint* position = (sint*) _mm_malloc(sizeof(int)*N, 64);
	int right_shift = (int)(log2(g_v_num))-8;
#pragma omp parallel num_threads(num_threads)
	{
		lint i, t; 
		sint threadid = omp_get_thread_num();
		lint p_per_thread = N/num_threads;
		if((p_per_thread * num_threads) != N) 
			p_per_thread++;
		lint start_n = p_per_thread*threadid;
		lint end_n = start_n + p_per_thread;
		if (end_n > N) 
			end_n = N;
		memset(histogram + threadid*NUMENTRIESPERTABLE, 0, sizeof(int)*NUMENTRIESPERTABLE);
		for (i = start_n; i < end_n; i++) {
			lint finalHash = idx_tmp[i] >> right_shift;
			position[i] = histogram[threadid*NUMENTRIESPERTABLE + finalHash];
			histogram[threadid*NUMENTRIESPERTABLE + finalHash]++;
		}
		//__BARRIER__;
		TREE_BARRIER(&tree_barrier, threadid, num_threads);
#pragma omp single 
		{
			int current_sum = 0;
			int prev_sum = 0;
			memset(bin_idx, 0, sizeof(int)*(NUMENTRIESPERTABLE+1));
			bin_idx[0] = 0;
			for (i = 0; i < NUMENTRIESPERTABLE; i++) {
				for (t = 0; t < num_threads; t++) {
					//cumsum[i] += histogram[t*NUMENTRIESPERTABLE + i-1];
					current_sum = prev_sum + histogram[t*NUMENTRIESPERTABLE + i];
					histogram[t*NUMENTRIESPERTABLE + i] = prev_sum;
					prev_sum = current_sum;
				}
				//cumsum[i] += cumsum[i-1];
				bin_idx[i+1] = current_sum;
			}
			if(current_sum != N) 
				ERROR_PRINT();
		}
		//__BARRIER__;
		TREE_BARRIER(&tree_barrier, threadid, num_threads);
		for (i = start_n; i < end_n; i++) {
			lint finalHash = idx_tmp[i] >> right_shift;
			lint newIndex = histogram[threadid*NUMENTRIESPERTABLE + finalHash] + position[i];
			//assign values
			tmp_diff_assign[newIndex] = diff_assign[i];
			tmp_pidx[newIndex] = pidx_tmp[i];
			tmp_idx[newIndex] = idx_tmp[i];
			tmp_node_idx[newIndex] = node_idx_tmp[i];
		}

	}
	//assign back
	lint i;
//#pragma omp parallel for num_threads(num_threads)
	for(i=0;i<size; i++){
		diff_assign[i] = tmp_diff_assign[i];
		pidx_tmp[i] = tmp_pidx[i];
		idx_tmp[i] = tmp_idx[i];
		node_idx_tmp[i] = tmp_node_idx[i];
	}
	_mm_free(histogram);
	_mm_free(position);
}

void to_msg(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank,
		t_buff *buf){
	//TQ_LOCKTYPE lock;
	//TQ_LOCKINIT(lock); 
	omp_lock_t writelock;
	omp_init_lock(&writelock);
	lint i,j=0;
	sint myrk=0;
	for(i = 0;i<num_nodes;i++){
		msg->idx_send[i] = 0;
	}
	lint lock_num=0;
	lint **msg_idx_remain = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);	
	lint **remain_idx = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);	
	for(i=0;i<num_threads;i++){
		msg_idx_remain[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		remain_idx[i] = (lint*)_mm_malloc(sizeof(lint)*num_nodes, 64);
		for(j=0;j<num_nodes;j++)
			msg_idx_remain[i][j]=0;
	}
	//prepare for the arrays that used in the synchronization
#ifdef USE_OMP
	omp_set_num_threads(num_threads);	
#endif
	//start computation
#ifdef USE_OMP
#pragma omp parallel shared(msg, buf, lock_num, msg_idx_remain, remain_idx) private(i, j, myrk)
	{
		myrk = omp_get_thread_num();
#endif
		if(myrk==0)
			tic_sr(t_list, 5);
		for(i=myrk;i<gs->v_size;i+=num_threads){
			for(j=(i==0?0:gs->vet_idx[i-1]); j<gs->vet_idx[i]; j++){
				lint source = i;
				lint target = gs->edge_idx[j];
				double target_val = diff[source];
				lint from = myrank;
				lint to = gs->edge_info[j].node_id;
				//if(to==0 && from ==3 && target == 1452)
				//	printf("ok the source is %lld in node %lld and the value is %lf\n", source, from, diff[source]);
				if(from!=to){
#ifdef USE_OMP
					buf->int_buffer[myrk][to][buf->buff_idx[myrk][to]]=target;
					buf->double_buffer[myrk][to][buf->buff_idx[myrk][to]]=target_val;
					buf->buff_idx[myrk][to] += 1;
					msg_idx_remain[myrk][to] += 1;
					if(buf->buff_idx[myrk][to]==buff_size){
						omp_set_lock(&writelock);
						//TQ_LOCK(lock);
						//printf("msg size is %lld buff size is %d|%d\n", msg->idx_send[to], buff_idx[myrk][to], buff_size);
						msg->idx_send[to]=cpy_from_buff_to_msg(msg, buf->int_buffer, buf->double_buffer, myrk, to, buf->buff_idx[myrk][to], msg->idx_send[to]);
						buf->buff_idx[myrk][to]=0;
						msg_idx_remain[myrk][to]=0;
						lock_num++;
						//TQ_UNLOCK(lock);
						omp_unset_lock(&writelock);
					}
#else
					lint idx =  msg->idx_send[to];
					msg->double_send[to][idx]=diff[source];
					msg->int_send[to][idx]=target;
					msg->idx_send[to]+=1;
#endif
				}
			}
		}		
		if(myrk==0)
			toc_sr(t_list, 5);
#ifdef USE_OMP
		//flush buffer, this part can be parallelized as well
#pragma omp barrier
		if(myrk==0)
			tic_sr(t_list, 6);
		//do some reduction work
		if(myrk==0){
			for(i=0;i<num_nodes;i++){
				lint sum = msg->idx_send[i];
				remain_idx[0][i]= msg->idx_send[i];
				for(j=1;j<num_threads;j++){
					sum += msg_idx_remain[j-1][i];
					remain_idx[j][i]= sum;
				}
			}
			//for(i=0;i<num_nodes;i++){
			//	for(j=0;j<num_threads;j++)
			//		printf("%lld|%lld|%lld ", remain_idx[j][i],  msg_idx_remain[j][i], buf->buff_idx[j][i]);
			//	printf("\n");
			//}
		}

#pragma omp barrier
		//move elements
		for(i=0;i<num_nodes;i++){
			lint to = i;
			omp_set_lock(&writelock);
			//printf("idx msg: %lld, idx buf: %lld\n", msg->idx_send[to], buf->buff_idx[myrk][to]);
			//if(myrank==0)
			//printf("idx msg: %lld, idx buf: %lld\n", remain_idx[myrk][to], msg_idx_remain[myrk][to]);
			//msg->idx_send[to]=cpy_from_buff_to_msg(msg, buf->int_buffer, buf->double_buffer, myrk, to, buf->buff_idx[myrk][to], msg->idx_send[to]);
			cpy_from_buff_to_msg(msg, buf->int_buffer, buf->double_buffer, myrk, to, msg_idx_remain[myrk][to], remain_idx[myrk][to]);
			buf->buff_idx[myrk][to]=0;
			omp_unset_lock(&writelock);
		}
		/*if(myrk==0){
			for(i=0;i<num_nodes;i++)
				msg->idx_send[i] = msg_idx_remain[num_threads-1][i]+remain_idx[num_threads-1][i];
		}*/
		/*char f_buf[100];
		sprintf(f_buf, "./data/log/logg_%d", myrank);
		FILE *f = fopen(f_buf, "w");
		if(myrk == 0){
			for(i=0;i<num_nodes;i++){
				for(j=0;j<num_threads;j++){
					lint to = i;
					//cpy_from_buff_to_msg(msg, buf->int_buffer, buf->double_buffer, j, to, msg_idx_remain[j][to], remain_idx[j][to]);
					sint id_buf = msg_idx_remain[j][to];
					lint id_msg = remain_idx[j][to];
					sint k;
					//if(myrank==0)
					//	printf("idx msg: %lld, idx buf: %lld\n", msg->idx_send[to], buf->buff_idx[j][to]);
					//msg->idx_send[to]=cpy_from_buff_to_msg(msg, buf->int_buffer, buf->double_buffer, j, to, buf->buff_idx[j][to], msg->idx_send[to]);
					//sint id_buf = buf->buff_idx[j][to];
					//lint id_msg = msg->idx_send[to];
					//if(myrank==0)
					fprintf(f, "idx msg: %lld, idx buf: %lld\n", id_msg, id_buf);
					for(k=0; k<id_buf; k++){
						msg->int_send[to][id_msg]=buf->int_buffer[j][to][k];
						msg->double_send[to][id_msg]=buf->double_buffer[j][to][k];
						id_msg++;
					}
					msg->idx_send[to]=id_msg;
					buf->buff_idx[j][to]=0;
				}
			}
		}
		fclose(f);*/
#pragma omp barrier
		if(myrk==0)
			toc_sr(t_list, 6);
	}
	//printf("num of locks used %lld \n",lock_num);
#endif
}

	lint 
cpy_from_buff_to_msg(t_msg *msg, lint ***int_buffer, 
		double ***double_buffer, sint myrk, 
		lint to, sint idx_buff, 
		lint idx_msg)
{
	lint i;
	sint id_buf = idx_buff;
	lint id_msg = idx_msg;
	for(i=0; i<id_buf; i++){
		msg->int_send[to][id_msg]=int_buffer[myrk][to][i];
		msg->double_send[to][id_msg]=double_buffer[myrk][to][i];
		id_msg++;
	}
	return id_msg;
}

void
reduce_msg(t_msg *msg, sint myrank, 
		lint **send_ub){
	int i,j;
	//there is a gather step
	//sort first
	for(i=0;i<num_nodes;i++){
		if(i==myrank)
			continue;
		q_sort_int(msg->int_send[i], msg->double_send[i], 0, send_ub[myrank][i]-1);
	}
}

void print_msg(double *diff, t_msg *msg, 
		t_csr *gs, sint myrank){
	lint i,j;
	printf("message copied: \n");
	for(i=0;i<num_nodes;i++)
	{
		for(j=0;j<msg->idx_send[i];j++)
			printf("%lld ", msg->int_send[i][j]);
		printf("\n");
	}	
}

void print_msg_file_send(FILE *f, double *diff, t_msg *msg, 
		t_csr *gs, sint myrank){
	lint i,j;
	fprintf(f, "==============message send===============: \n");
	for(i=0;i<num_nodes;i++)
	{
		fprintf(f, "num msg %lld ", msg->idx_send[i]);
		for(j=0;j<msg->idx_send[i];j++)
			fprintf(f, "%lld|%lf ", msg->int_send[i][j], msg->double_send[i][j]);
		fprintf(f, "\n");
	}	
}
void print_msg_file_recv(FILE *f, double *diff, t_msg *msg, 
		t_csr *gs, sint myrank){
	lint i,j;
	//#ifdef USE_DEBUG
	fprintf(f, "==============message recv===============: \n");
	for(i=0;i<num_nodes;i++)
	{
		for(j=0;j<msg->idx_recv[i];j++)
			fprintf(f, "%lld|%lf ", msg->int_recv[i][j], msg->double_recv[i][j]);
		fprintf(f, "\n");
	}	
	//#endif
}

void print_sendub_file(FILE *f, lint**send_ub){
	lint i,j;
	for(i=0;i<num_nodes;i++){
		for(j=0;j<num_nodes;j++)
			DPRINTF(3, "%lld ", send_ub[i][j]);
		DPRINTF(3, "\n");
	}
}

void pp_msg(t_msg *msg, t_csr *gr, 
		sint myrank, lint **send_ub, double purpose_jump)
{
	lint i,j;
	for(i=0;i<num_nodes;i++){
		for(j=0;j<msg->idx_recv[i];j++){
			lint idx = msg->int_recv[i][j] - gr->offset + OFFSET;
			if(myrank==0 && idx == 0)
				printf("%lld %lf\n", i, msg->double_recv[i][j]);
		}
	}	
}

void from_msg(t_msg *msg, t_csr *gr, 
		sint myrank, lint **send_ub, 
		double purpose_jump, lint **idx_start,
		lint **idx_end)
{
	lint i,j;
	//partition work for every thread
	for(i=0;i<num_nodes;i++)
		for(j=0;j<num_threads;j++){
			idx_start[i][j]=0;
			idx_end[i][j]=0;
		}
	for(i=0;i<num_nodes;i++){
		if(i==myrank)
			continue;
		lint avg = msg->idx_recv[i]/num_threads+1;
		//printf("avg %lld\n", avg);
		for(j=0;j<num_threads;j++){
			idx_start[i][j]=j*avg;
			idx_end[i][j]=(j==num_threads-1?msg->idx_recv[i]:(j+1)*avg);
			if(j>0){
				lint k=j*avg-1;
				while(1){
					if(msg->int_recv[i][k]==msg->int_recv[i][k+1]){
						k--;
						idx_start[i][j]--;
						idx_end[i][j-1]--;
					}
					else
						break;
				}
			}
		}
	}
	//char f_buf[100];
	//sprintf(f_buf, "./data/log/log_%d", myrank);
	//FILE *logf = fopen(f_buf, "w"); 
	//for(i=0;i<num_nodes;i++){
	//	for(j=0;j<num_threads;j++)
	//		fprintf(logf, "%lld|%lld ", idx_start[i][j], idx_end[i][j]);
	//	fprintf(logf, "\n");
	//}
	//for(i=0;i<msg->idx_recv[1];i++)
	//	fprintf(logf, "%lld %lld %lf\n", i, msg->int_recv[1][i], msg->double_recv[1][i]);
	//fclose(logf);
	//exit(1);
	//	lint k;
	//	for(i=0; i<num_nodes; i++)
	//		for(j=0;j<num_threads;j++)
	//			for(k=0;k<gr->v_size; k++)
	//				dirty[i][j][k]=FALSE;

#ifdef USE_OMP
	omp_set_num_threads(num_threads);
#pragma omp parallel shared(msg, gr, purpose_jump) private(i,j) 
	{
		sint thread_id = omp_get_thread_num();
#else
		sint thread_id=0;
#endif

		for(i=0;i<num_nodes;i++){
			lint start = idx_start[i][thread_id];
			lint end = idx_end[i][thread_id];
			for(j=start;j<end;j++){
				lint idx = msg->int_recv[i][j] - gr->offset + OFFSET;
				gr->vet_info[idx].wb += purpose_jump*msg->double_recv[i][j];
				//dirty[i][thread_id][idx]=TRUE;
				//if(myrank==3 && idx ==40856)
				//	printf("I am %d, I am updating %lld|%lld from node %lld the value is %lf\n",  myrank, msg->int_recv[i][j],idx, i, msg->double_recv[i][j]);
			}
#pragma omp barrier
		}	
#ifdef USE_OMP
	}
#endif
	/*for(i=0;i<num_nodes;i++){
	  if(i==myrank)
	  continue;
	  for(k=0;k<gr->v_size;k++){
	  sint count=0;
	  for(j=0;j<num_threads;j++){
	  if(dirty[i][j][k]==TRUE)
	  count++;
	  }
	  if(count>1)
	  printf("hey here is wrong!");
	  }
	  }*/
	/*for(i=0;i<num_nodes;i++){
	  if(i==myrank)
	  continue;
	  lint start = 0;
	  lint end = msg->idx_recv[i];
	  for(j=start;j<end;j++){
	  lint idx = msg->int_recv[i][j] - gr->offset + OFFSET;
	  gr->vet_info[idx].wb += purpose_jump*msg->double_recv[i][j];
	  if(myrank==0 && idx==1)
	  printf("%lf \n", msg->double_recv[i][j]);
	//if(myrank==9 && idx == 47245)
	//	printf("I am %d, I am updating %lld|%lld from node %lld the value is %lf\n",  myrank, msg->int_recv[i][j], idx, i, msg->double_recv[i][j]);
	}
	}*/
}

void run_pagerank_csr_mpi(char *file, char *dir)
{
	t_list = (t_time_list*)malloc(sizeof(t_time_list));
	t_list->list = (t_time_elem*)malloc(sizeof(t_time_elem)*10);
	init_timer(t_list, 10);
	lint i=0,j=0;
	t_buff *buf = (t_buff*)malloc(sizeof(t_buff));
	buf->size=buff_size;
	buf->buff_idx = (sint**)malloc(sizeof(sint*)*num_threads);
	buf->double_buffer = (double***)malloc(sizeof(double**)*num_threads);
	buf->int_buffer = (lint***)malloc(sizeof(lint**)*num_threads);
	init_t_buff(buf, buf->size);
	//data structure for from_msg
	lint **idx_start = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	lint **idx_end = (lint**)_mm_malloc(sizeof(lint*)*num_nodes, 64);
	for(i=0;i<num_nodes;i++){
		idx_start[i] = (lint*)_mm_malloc(sizeof(lint)*num_threads, 64);
		idx_end[i] = (lint*)_mm_malloc(sizeof(lint)*num_threads, 64);
	}

	//double rand_jump = RDM_JMP/(u/t_msouble)graph->v_size;
	double rand_jump = RDM_JMP;
	double purpose_jump = 1-RDM_JMP;
	lint iter=0;

	/* Initialize mpi */
	int *argc = (int*)malloc(sizeof(int));
	char ***argv = (char***)malloc(sizeof(char**));
	//MPI_Init(argc, argv);	
	sint myrank, size;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	// init files
	char f_buf[100];

	sprintf(f_buf, "%slog/%s_log_%d_%d",dir, file, myrank, num_nodes);
	FILE *lg = fopen(f_buf, "w");
	//if((lg = fopen(f_buf, "w"))==NULL){
	//	printf("the file %s does not exists!\n", f_buf);
	//	ERROR_PRINT();
	//}
	////////////////////////////////
	lint ii, jj;
	lint upper_bound, lower_bound;
	/* reading graph, init some parameters */
	t_csr *gs = (t_csr*)malloc(sizeof(t_csr));
	sprintf(f_buf, "%s%s_sd_%d_%d", data_extend, file, myrank, num_nodes);
	if(bin==FALSE){
		scan_csr_idx(gs, f_buf, dir, SEND);
		read_graph_csr(gs, f_buf, dir, SEND);
	}
	else{
		read_csr_bin(gs, f_buf,dir, SEND, myrank);
	}
	t_csr *gr = (t_csr*)malloc(sizeof(t_csr));
	sprintf(f_buf, "%s%s_rc_%d_%d", data_extend, file, myrank, num_nodes);
	if(bin==FALSE){
		scan_csr_idx(gr, f_buf, dir, RECV);
		read_graph_csr(gr, f_buf, dir, RECV);
	}
	else{
		read_csr_bin(gr, f_buf,dir, RECV, myrank);
	}
	lint *distrib = (lint*)malloc(sizeof(lint)*num_nodes);
	if(bin==FALSE){
		read_recip(file, dir, distrib, gr, myrank);
		if(myrank==0)
			printf("%1.10f\n", gr->vet_info[0].recip);
	}
	else{
		read_recip_bin(file, dir, distrib, gr, myrank);
	}
	init_node_id(gs, distrib, num_nodes);
	init_weight(gr);
	t_msg *msg = (t_msg*)malloc(sizeof(t_msg));
	lower_bound = gr->offset-OFFSET;
	upper_bound = gr->v_size+gr->offset-OFFSET;
	//printf("%d %d %d \n", myrank, lower_bound, upper_bound);
	double *diff_pr = (double*)_mm_malloc(sizeof(double)*gr->v_size, 64);
	double *diff_tmp = (double*)_mm_malloc(sizeof(double)*gs->e_size, 64);
	lint *pidx_tmp = (lint*)malloc(sizeof(lint)*gs->e_size);
	lint *diff_assign = (lint*)malloc(sizeof(lint)*gs->e_size);
	lint *idx_tmp = (lint*)malloc(sizeof(lint)*gs->e_size);
	lint *node_idx_tmp = (lint*)malloc(sizeof(lint)*gs->e_size);
	lint *bin_idx = (lint*)malloc(sizeof(lint)*(1+HASH_BINS));
	for(i=0;i<gr->v_size;i++){
		diff_pr[i]=gr->vet_info[i].weight * gr->vet_info[i].recip;
	}

	sprintf(f_buf, "%sparex/%s_pidx_%d_%d", dir, file, myrank, num_nodes);
	FILE *reader_pidx;
	if((reader_pidx=fopen(f_buf, "r")) == NULL){
		printf("the file %s does not exists!\n", f_buf);
		ERROR_PRINT();
	}
	for(i =0; i< gs->e_size; i++){
		lint pval;
		fscanf(reader_pidx, "%lld\n", &pval);
		pidx_tmp[i] = pval;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(j=0; j<gs->e_size; j++){
		diff_assign[j] = j;
		idx_tmp[j] = gs->edge_idx[j];
		node_idx_tmp[j] = gs->edge_info[j].node_id;
	}
	//for(i=0;i<gs->e_size; i++)
	//	printf("%lld %lld\n", idx_tmp[i], pidx_tmp[i]);
	//read pidx_tmp
	//partition the msg
	partition_msg(diff_assign, pidx_tmp, idx_tmp, bin_idx, node_idx_tmp, gs->e_size);
	reverse_diff_assign(diff_assign, gs->e_size);
	//sprintf(f_buf, "data/rand/log/2_20_l_%d", myrank);
	//FILE *l = fopen(f_buf, "w");
	//for(i=0;i<gs->e_size; i++)
	//	fprintf(l, "%lld\n", diff_assign[i]);
	//fprintf(l, "\n");
	//fclose(l);
	//printf("finished partition\n");
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	lint **send_ub = (lint**)malloc(sizeof(lint*)*num_nodes);
	for(i=0;i<num_nodes;i++)
		send_ub[i] = (lint*)malloc(sizeof(lint)*num_nodes);
	if(bin==FALSE)
		read_send_ub(file, dir, send_ub);
	else
		read_send_ub_bin(file, dir, send_ub);
	//if(myrank==0){
	//	printf("the send ub is\n");
	//	for(i=0;i<num_nodes; i++){
	//		for(j=0;j<num_nodes; j++)
	//			printf("%lld ", send_ub[i][j]);
	//		printf("\n");
	//	}
	//}


	init_msg(msg, send_ub, myrank);
	for(ii=0;ii<gr->v_size;ii++)
		gr->vet_info[ii].wb = rand_jump;

	/* major computation */
#ifdef USE_OMP
	omp_set_num_threads(num_threads);	
#endif
	sint terminate = FALSE;
	sint term_buf[10];
	sint term[num_nodes];
	while(terminate  == FALSE){
		for(ii=0;ii<num_nodes;ii++)
			term[ii]=TRUE;
		iter++;
		terminate = TRUE;
		if(myrank==0){
			DPRINTF(1, "ITER %lld ", iter);
			for(ii=0;ii<10;ii++){
				DPRINTF(1, "%f|%f ", gr->vet_info[ii].wb, gr->vet_info[ii].recip);
			}
			//DPRINTF(1, "%f|%f ", gr->vet_info[47245].wb, gr->vet_info[53765].recip);
			DPRINTF(1, "\n");
		}
		tic_sr(t_list, 0);
#ifdef USE_OMP
#pragma omp parallel for shared (gs, gr)
#endif
		for(ii=0;ii<gr->v_size;ii++){
			for(jj=(ii==0?0:gr->vet_idx[ii-1]); jj<gr->vet_idx[ii]; jj++){
				lint source = ii;
				lint target = gr->edge_idx[jj];
				if(target>= lower_bound && target<upper_bound){
					gr->vet_info[source].wb+=purpose_jump*diff_pr[target - gr->offset + OFFSET];
					// prefetch next target value
					lint target_next = gr->edge_idx[jj+prefetch_dis];
					_mm_prefetch((char *)(diff_pr + target_next), _MM_HINT_T0);
				}
			}
		}
		toc_sr(t_list, 0);
		for(ii=0;ii<num_nodes;ii++)
			msg->idx_recv[ii]=0;
		DPRINTF(3, "node %d, on iter %lld, finished local computation!\n", myrank, iter);
		//communication send 
		sint num_comm = log(num_nodes)/log(2);
		sint ic, it;
		sint bucket[num_nodes];
		sint role[num_nodes];
		tic_sr(t_list, 1);
#ifdef USE_SORTMSG	
		to_msg_no_buf(diff_pr, msg, gs, myrank);
		reduce_msg(msg, myrank, send_ub);
#else
		assign_diff(diff_tmp, diff_assign, diff_pr, gs, myrank);
		set_msg(msg, send_ub, myrank);
		to_msg_reduce(diff_tmp, msg, pidx_tmp, idx_tmp, bin_idx, node_idx_tmp);
#endif
		toc_sr(t_list, 1);

		MPI_Request request_send_int[num_nodes];
		MPI_Request request_send_double[num_nodes];
		MPI_Request request_recv_int[num_nodes];
		MPI_Request request_recv_double[num_nodes];
		MPI_Status status;

		tic_sr(t_list, 2);
		sint msg_iter=1;
#ifdef USE_SYNC
		for(it=0;it<msg_iter;it++){
			for(ic=0; ic<num_comm; ic++){
				compute_bucket_role(bucket, role, ic);
				if(role[myrank]==SEND_FIRST){
					for(ii=0; ii<num_nodes; ii++){
						if(bucket[ii]!=bucket[myrank] || role[ii]!=RECV_FIRST){
							continue;
						}
						//odd send to even then recv from even
						lint send_avg = send_ub[myrank][ii]/msg_iter+1;
						lint send_start = send_avg*it;
						lint send_size = (send_ub[myrank][ii]-send_start)>=send_avg?send_avg:(send_ub[myrank][ii]-send_start);
						//printf("myrank %d, to %lld msg_iter %d avg:%lld start:%lld size:%lld ub: %lld\n", myrank, ii, it, send_avg, send_start, send_size, send_ub[myrank][ii]);
						MPI_Send(msg->double_send[ii]+send_start, send_size, MPI_DOUBLE,  ii, 2, MPI_COMM_WORLD );
						MPI_Send(msg->int_send[ii]+send_start, send_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD);
						lint recv_avg = send_ub[ii][myrank]/msg_iter+1;
						lint recv_start = recv_avg*it;
						lint recv_size = (send_ub[ii][myrank]-recv_start)>=recv_avg?recv_avg:(send_ub[ii][myrank]-recv_start);
						MPI_Recv(msg->double_recv[ii]+recv_start, recv_size, MPI_DOUBLE, ii, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Recv(msg->int_recv[ii]+recv_start, recv_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						msg->idx_recv[ii]=send_ub[ii][myrank];
						//if(iter==1 && myrank==0 && ii==1)
						//	printf("send ub %lld %lld\n", send_ub[ii][myrank], recv_size);
						//if(myrank==0 && ii==1 && iter ==1)
						//	for(i=0;i<recv_size;i++)
						//		printf("%d %lld\n", i+recv_start, msg->int_recv[ii][i+recv_start]);
					}
				}
				else if(role[myrank]==RECV_FIRST){
					for(ii=0; ii<num_nodes; ii++){
						if(bucket[ii]!=bucket[myrank] || role[ii]!=SEND_FIRST){
							continue;
						}
						//even recv from odd then send to odd
						lint recv_avg = send_ub[ii][myrank]/msg_iter+1;
						lint recv_start = recv_avg*it;
						lint recv_size = (send_ub[ii][myrank]-recv_start)>=recv_avg?recv_avg:(send_ub[ii][myrank]-recv_start);
						//if(myrank==1 && ii==0)
						//printf("myrank %d, to %lld msg_iter %d avg:%lld start:%lld size:%lld ub: %lld\n", myrank, ii, it, recv_avg, recv_start, recv_size, send_ub[ii][myrank]);
						MPI_Recv(msg->double_recv[ii]+recv_start, recv_size, MPI_DOUBLE, ii, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Recv(msg->int_recv[ii]+recv_start, recv_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						lint send_avg = send_ub[myrank][ii]/msg_iter+1;
						lint send_start = send_avg*it;
						lint send_size = (send_ub[myrank][ii]-send_start)>=send_avg?send_avg:(send_ub[myrank][ii]-send_start);
						MPI_Send(msg->double_send[ii]+send_start, send_size, MPI_DOUBLE,  ii, 2, MPI_COMM_WORLD);
						MPI_Send(msg->int_send[ii]+send_start, send_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD);
						msg->idx_recv[ii]=send_ub[ii][myrank];
					}
				}
				//MPI_Barrier(MPI_COMM_WORLD);
			}
		}
#else
		for(it=0;it<msg_iter;it++){
			for(ii=0;ii<num_nodes;ii++){
				if(ii==myrank)
					continue;
				lint send_avg = send_ub[myrank][ii]/msg_iter+1;
				lint send_start = send_avg*it;
				lint send_size = (send_ub[myrank][ii]-send_start)>=send_avg?send_avg:(send_ub[myrank][ii]-send_start);
				MPI_Isend(msg->double_send[ii]+send_start, send_size, MPI_DOUBLE,  ii, 2, MPI_COMM_WORLD, &(request_send_double[ii]));
				MPI_Isend(msg->int_send[ii]+send_start, send_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD, &(request_send_int[ii]));
			}
		}
		for(it=0;it<msg_iter;it++){
			for(ii=0;ii<num_nodes;ii++){
				if(ii==myrank)
					continue;
				lint recv_avg = send_ub[ii][myrank]/msg_iter+1;
				lint recv_start = recv_avg*it;
				lint recv_size = (send_ub[ii][myrank]-recv_start)>=recv_avg?recv_avg:(send_ub[ii][myrank]-recv_start);
				MPI_Irecv(msg->double_recv[ii]+recv_start, recv_size, MPI_DOUBLE, ii, 2, MPI_COMM_WORLD, &(request_recv_double[ii]));
				MPI_Irecv(msg->int_recv[ii]+recv_start, recv_size, MPI_LONG_LONG, ii, 0, MPI_COMM_WORLD, &(request_recv_int[ii]));
				msg->idx_recv[ii]=send_ub[ii][myrank];
			}
			//printf("come here\n");
		}
		//printf("before log\n");
		for(ii=0;ii<num_nodes;ii++){
			if(ii == myrank)
				continue;
			MPI_Wait(&(request_send_double[ii]), &status);
			MPI_Wait(&(request_recv_double[ii]), &status);
			MPI_Wait(&(request_send_int[ii]), &status);
			MPI_Wait(&(request_recv_int[ii]), &status);
		}
		//printf("see log!\n");
		//for(ii=0;ii<num_nodes;ii++){
		
		//}
#endif
		//if(iter==1){
		//	sprintf(f_buf, "/home/yinzhaom/SNY/data/log/log_%d", myrank);
		//	FILE *logfile = fopen(f_buf, "w");
		//	int node, vet;
		//	for(node=0;node<num_nodes;node++){
		//		if(node==myrank)
		//			continue;
		//		fprintf(logfile, "from node %d\n", node);
		//		for(vet=0;vet<send_ub[node][myrank];vet++)
		//			fprintf(logfile, "%d %lld %lf\n", vet, msg->int_recv[node][vet], msg->double_recv[node][vet]);
		//	}
		//}
		//printf("finished msg passing of rank %d\n", myrank);
		toc_sr(t_list, 2);
		if(myrank==0){
			for (jj=0;jj<num_nodes;jj++)
				DPRINTF(3, "%lld ", msg->idx_recv[jj]);
			DPRINTF(3, "\n");
		}
		/////////////communication recv////////////////////
		tic_sr(t_list, 3);
		from_msg(msg, gr, myrank,send_ub, purpose_jump, idx_start, idx_end);
		toc_sr(t_list, 3);
		//printf("FINIHSED from msg of rank %d\n", myrank);
		//if(iter==1)
		//	pp_msg(msg, gr, myrank,send_ub, purpose_jump);
		//////////////////check termination//////////////
		tic_sr(t_list, 4);
		for(ii=0;ii<gr->v_size;ii++){
			diff_pr[ii]=(gr->vet_info[ii].wb - gr->vet_info[ii].weight)*gr->vet_info[ii].recip;
		}
		for(ii=0;ii<gr->v_size;ii++){
			if(fabs(diff_pr[ii])>THRESH)
				terminate = FALSE;
			gr->vet_info[ii].weight = gr->vet_info[ii].wb;
		}
		term_buf[0]=terminate;
		for(ii=0; ii<num_nodes; ii++){
			if(ii==myrank)
				continue;
			MPI_Send(term_buf, 10, MPI_INT, ii, 1, MPI_COMM_WORLD);
		}
		for(ii=0; ii<num_nodes; ii++){
			if(ii==myrank)
				continue;
			MPI_Recv(term_buf, 10, MPI_INT, ii, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			term[ii]=term_buf[0];
		}
		term[myrank] = terminate;
		for(ii=0;ii<num_nodes;ii++)
			if(term[ii]==FALSE)
				terminate = FALSE;
		if(terminate == TRUE)
			break;	
		toc_sr(t_list, 4);
		//printf("finished check termination of rank %d\n", myrank);
	}
	if(myrank==MASTER){
		DPRINTF(1, "Final pagerank values: \n");
		for(ii=0;ii<10;ii++){
			DPRINTF(1, "%lld|%lf ", (ii+gr->offset), gr->vet_info[ii].weight);
		}
		DPRINTF(1, "\n");
	}
	if(myrank==MASTER){
		printf("iterations taken: %lld\n", iter);
		lint recv_size=0;
		for(i=0;i<num_nodes;i++){
			if(i==myrank)
				continue;
			recv_size += send_ub[i][myrank];
		}
		print_result(t_list,  iter, gs->e_size, recv_size);
	}
	/* Shut down mpi */
	MPI_Finalize();
	//fclose(lg);
} 

void 
compute_bucket_role(sint *bucket, sint *role, sint ii){
	sint i,j=0;
	sint buck_id=0;
	sint num_bucket = (sint)pow(2, ii);
	sint bucket_size = (sint)num_nodes/num_bucket;
	sint half_size = (sint)bucket_size/2;
	sint low=0;
	for(i=0;i<num_bucket;i++){
		for(j=0;j<bucket_size;j++){
			sint id = i*bucket_size + j;
			bucket[id] = i;
			if(j<half_size)
				role[id]=SEND_FIRST;
			else
				role[id]=RECV_FIRST;
		}
	}
}
