#include "graph.h"
#include "array.h"
#include "list.h"
#include "timer.h"

t_time_list *t_list_sssp;
lint **b_start;
lint **b_size;
lint *r_start;
lint *r_size;
lint *s_start;
lint *s_size;
lint **local_rmv_cnt;
lint **local_add_cnt;
lint iter = 0;
t_csr *gl;
t_csr *gh;
lint *bin_idx;
sint is_heavy;
p_arr tmp_B;
lint ***rev_stats;
int inner_i=0;
#define SUC 0
#define CFLT 1
#define BACK 2
#define SAME 3
#define TTL 4
//int prefetch_dis = 1;
//lint *local_b_total;

//this is for serial code
void reorder_arr(p_arr arr, t_csr *gs, double *tmp);
void			
relax_arr_by_partition(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
relax_arr_single(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
request_arr_single(t_csr *gs, p_arr arr, p_arr R);
void
remember_arr_single(p_arr B, p_arr S, 
		lint *b_total, p_vidx vidx);

//this is for parallel code
void
print_sssp_results(p_time_list pl);
void 
compute_request_partition(p_arr R);
void 
print_B_map(FILE *log, p_arr *B, 
		lint iter);
void
set_vidx(p_arr B, p_vidx vidx, sint iter);
int
check_B_total_on_iter(p_arr **B, lint iter);
void 
compute_request_idx(sint thread_id, t_csr *gs, t_arr *arr, FILE *log);
void 
reduce_request_idx();
void
balance_request(p_arr RR, p_arr R, 
		sint thread_id);
lint
find_edge_idx_arr(t_csr *gs, lint source, lint target);
void
relax_arr(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx);
void
request_arr(t_csr *gs, p_arr arr, p_arr R);
void
remember_arr(p_arr B, p_arr S, lint *b_total, p_vidx vidx);
lint 
check_total(lint *b_total);
void recover_arr(p_arr *B, p_vidx vidx, 
		lint iter);
void
compute_RR_partition(p_arr RR);

//void
//compute_balance_index(p_arr *B, lint **balance_idx_start, 
//			lint**balance_idx_end){
//	int i, j;
//	lint total=0;
//	lint avg=0;
//	sint num_less=0;
//	lint total_less=0;
//	double less_weight[num_threads];
//	//compute how much work should be assigned to which thread
//	for(i=0;i<num_threads;i++){
//		total+=B[i]->num;
//		less_weight[i]=0;
//	}	
//	avg = total/num_threads;
//	for(i=0;i<num_threads;i++){
//		if(B[i]->num < avg)
//			total_less += avg-B[i]->num;
//	}
//	for(i=0;i<num_threads;i++){
//		if(B[i]->num < avg)
//			less_weight[i] =  (double)(avg - B[i]->num)/(double)total_less;
//	}
//	//compute the exact index
//	for(i=0;i<num_threads;i++){	
//		lint off = avg;
//		lint more = B[i]->idx - avg;
//		if(more > 0){
//			for(j=0; j<num_threads; j++){
//				if(B[j]->idx < avg){
//					balance_idx_start[i][j] = off;
//					lint deficit = more*less_weight[j];
//					off += deficit;
//					balance_idx_end[i][j] = off;
//					B[i]->idx -= deficit; 
//				}
//			}
//		}
//	}
//}

//void 
//balance(p_arr *B, lint **balance_idx_start,
//		lint**balance_idx_end, sint thread_id){	
//	int from, from_idx;
//	sint to = thread_id;
//	for(from=0; from<num_threads; from++){
//		for(from_idx=balance_idx_start[from][to]; from_idx < balance_idx_end[from][to]; from_idx++){
//			lint to_idx = B[to]->idx;
//			B[to]->arr[to_idx].int_val = B[from]->arr[from_idx].int_val;	
//			B[to]->arr[to_idx].int_val1 = B[from]->arr[from_idx].int_val1;	
//			B[to]->arr[to_idx].double_val = B[from]->arr[from_idx].double_val;	
//			B[from]->arr[from_idx].int_val = -1;
//			B[from]->arr[from_idx].int_val1 = -1;
//			B[from]->arr[from_idx].double_val = -1;
//			B[from]->num--;
//			B[to]->idx++;
//			B[to]->num++;
//		}
//	}
//}


void
run_sspath_arr(t_csr *gs, p_list *predecesor, 
		lint source, double delta){
	lint i, j;		
#ifdef REV_STATS
	rev_stats = (lint***)_mm_malloc(sizeof(lint**)*num_threads, 64);
	for(j=0; j<num_threads; j++){
		rev_stats[j] = (lint**)_mm_malloc(sizeof(lint*)*buck_size, 64);
		for(i=0; i<buck_size; i++)
			rev_stats[j][i] = (lint*)_mm_malloc(sizeof(lint)*5, 64);
	}
#endif
	char f_buf[100];
	lint sorting_thresh = 1024;
	sprintf(f_buf, "./data/log/plog");
	FILE *log = fopen(f_buf, "w");
	t_list_sssp = (t_time_list*)malloc(sizeof(t_time_list));
	t_list_sssp->list = (t_time_elem*)malloc(sizeof(t_time_elem)*10);
	init_timer(t_list_sssp, 10);
	
	//separate graph into light and heavy graphs
	gl = (t_csr*)malloc(sizeof(t_csr));
	gh = (t_csr*)malloc(sizeof(t_csr));
	g_v_num = gs->v_size;
	separate_csr(gs, gl, gh, delta);	
	//{
	//	lint start = gh->vet_idx[743786];
	//	lint end = gh->vet_idx[743787];
	//	printf("start %lld end %lld\n", start, end);
	//	for(i=start;i<end;i++)
	//		printf("target %lld\n", gh->edge_idx[i]);
	//	exit(1);
	//}
	//prepare for data structure
	lint arr_size = gs->v_size;
	//buck_size is how many buckets needed
	p_arr *B = (p_arr*)_mm_malloc(sizeof(p_arr*)*buck_size, 64);		
	for(i=0;i<buck_size;i++){
		B[i] = (p_arr)_mm_malloc(sizeof(t_arr), 64);
		create_arr(B[i], arr_size*2);
	}
	printf("finished preparing for the bucket!\n");
	tmp_B = (p_arr)_mm_malloc(sizeof(t_arr), 64);
	create_arr(tmp_B, arr_size);
	p_arr R = (p_arr)_mm_malloc(sizeof(t_arr), 64);
	double *tmp = (double*)_mm_malloc(sizeof(double)*gs->e_size, 64);
	p_arr tmp_R = (p_arr)_mm_malloc(sizeof(t_arr), 64);
	p_arr S = (p_arr)_mm_malloc(sizeof(t_arr), 64);
	p_vidx vidx = (p_vidx)_mm_malloc(sizeof(t_vidx)*arr_size, 64);
	for(i=0;i<arr_size;i++){
		vidx[i].buck_id=-1;
		vidx[i].buck_id=-1;
	}
	create_arr(R, gs->e_size);
	create_arr(tmp_R, gs->e_size);
	create_arr(S, gs->e_size);
	for(i=0;i<gs->v_size;i++)
		gs->vet_info[i].weight = INFINITY;
	//for transfering data to R
	b_start = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);
	b_size = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64);
	local_rmv_cnt = (lint**)_mm_malloc(sizeof(lint)*num_threads, 64);
	local_add_cnt = (lint**)_mm_malloc(sizeof(lint)*num_threads, 64);
	//local_b_total = (lint*)_mm_malloc(sizeof(lint)*num_threads, 64);
	for(i=0;i<num_threads;i++){
		b_start[i]=(lint*)_mm_malloc(sizeof(lint)*buck_size, 64);
		b_size[i]=(lint*)_mm_malloc(sizeof(lint)*buck_size, 64);
		local_rmv_cnt[i] = (lint*)_mm_malloc(sizeof(lint)*buck_size, 64);
		local_add_cnt[i] = (lint*)_mm_malloc(sizeof(lint)*buck_size, 64);
	}
	//for transfering data from R
	r_start = (lint*)_mm_malloc(sizeof(lint)*num_threads*64, 64);
	r_size = (lint*)_mm_malloc(sizeof(lint)*num_threads*64, 64);
	s_start = (lint*)_mm_malloc(sizeof(lint)*num_threads*64, 64);
	s_size = (lint*)_mm_malloc(sizeof(lint)*num_threads*64, 64);
	//relax_arr source vertex
	add_arr_two(B[0], source, 0.0);
	vidx[source].buck_id = 0;
	vidx[source].buck_pos = 0;
	gs->vet_info[source].weight = 0.0;
	lint b_total=1;

	bin_idx = (lint *)malloc(sizeof(lint)*(1+HASH_BINS));
	printf("finished preparation of all resources!\n");
	//major computing
#ifdef USE_OMP
	omp_set_num_threads(num_threads); 
#endif
	while(b_total>0){
		empty_arr(R);
		inner_i=0;
		while(B[iter]->num>0){
			//for(i=0;i<B[iter]->idx; i++)
			//	printf("requst %lld\n", B[iter]->arr[i].int_val);
			inner_i ++;
#ifdef USE_OMP
			tic_sr(t_list_sssp, 0);
			request_arr(gl, B[iter], tmp_R);
			toc_op(t_list_sssp, 0, B[iter]->num, R->num, 88, 64);
			//88 = 24 + 64
			//64 = 2*24 + 16
			//q_sort_arr(R, 0, R->num-1);
			tic_sr(t_list_sssp, 6);
			if(use_partition == TRUE && tmp_R->idx > sorting_thresh)
				partition(tmp_R, R, bin_idx);
			else{
				merge_sort(R, tmp_R);
				reorder_arr(R, gs, tmp);
				compute_request_partition(R);
			}
			toc_sr(t_list_sssp, 6);

			tic_sr(t_list_sssp, 1);
			lint os_num = S->num;
			remember_arr(B[iter], S, &b_total, vidx);
			toc_op(t_list_sssp, 1, S->num-os_num, 0, 136, 0);
			//136 = 24 + 2*24 + 64

			tic_sr(t_list_sssp, 2);
			//printf("iter %lld relax light\n", iter);
			is_heavy = FALSE;
			if(use_partition == TRUE && tmp_R->idx > sorting_thresh){
				relax_arr_by_partition(gs, predecesor, B, R, delta, &b_total, vidx);
			}
			else
				relax_arr(gs, predecesor, B, R, delta, &b_total, vidx);
			//printf("iter %lld:%d %p\n", iter, inner_i, B[4]);
		//	printf("iter %lld after relaxing light, there are:\n", iter);
		//	for(i=0;i<buck_size;i++){
		//		if(B[i]->num != 0)
		//			printf("B[%lld] has %lld elems\n", i, B[i]->num);
		//	}
			toc_op(t_list_sssp, 2, tmp_R->num, 0, 88, 0);
#else

			tic_sr(t_list_sssp, 0);
			request_arr_single(gl, B[iter], R);
			toc_op(t_list_sssp, 0, B[iter]->num, R->num, 88, 64);
			//88 = 24 + 64
			//64 = 2*24 + 16
			tic_sr(t_list_sssp, 6);
			//q_sort_arr(R, 0, R->num-1);
			//reorder_arr(R, gs);
			toc_sr(t_list_sssp, 6);

			tic_sr(t_list_sssp, 1);
			remember_arr_single(B[iter], S, &b_total, vidx);
			toc_op(t_list_sssp, 1, S->num, 0, 136, 0);
			//136 = 24 + 2*24 + 64
			//printf("726471-1 idx is %lld %lld\n", vidx[726471].buck_id, vidx[726471].buck_pos);
			//printf("b_total %lld\n", b_total);
			tic_sr(t_list_sssp, 2);
			relax_arr_single(gs, predecesor, B, R, delta, &b_total, vidx);
			toc_op(t_list_sssp, 2, R->num, 0, 88, 0);
			//printf("726471-2 idx is %lld %lld\n", vidx[726471].buck_id, vidx[726471].buck_pos);
			//printf("b_total %lld\n", b_total);
#endif
		}
		//printf("heavy iter +++++++++++++\n");
#ifdef USE_OMP
		tic_sr(t_list_sssp, 3);
		request_arr(gh, S, tmp_R);
		toc_op(t_list_sssp, 3, S->num, R->num, 88, 64);

		//q_sort_arr(R, 0, R->num-1);
		tic_sr(t_list_sssp, 6);
		if(use_partition == TRUE && tmp_R->idx > sorting_thresh)
			partition(tmp_R, R, bin_idx);
		else{
			merge_sort(R, tmp_R);
			reorder_arr(R, gs, tmp);
			compute_request_partition(R);
		}
		//partition(R, tmp_R, bin_idx);
		toc_sr(t_list_sssp, 6);

		tic_sr(t_list_sssp, 4);
		//printf("iter %lld relax heavy\n", iter);
		is_heavy = TRUE;
		if(use_partition == TRUE && tmp_R->idx > sorting_thresh)
			relax_arr_by_partition(gs, predecesor, B, R, delta, &b_total, vidx);
		else
			relax_arr(gs, predecesor, B, R, delta, &b_total, vidx);
		//printf("heavy iter %lld %p\n", iter, B[4]);
		//printf("iter %lld after relaxing heavy, there are:\n", iter);
		//for(i=0;i<buck_size;i++){
		//	if(B[i]->num != 0)
		//		printf("B[%lld] has %lld elems\n", i, B[i]->num);
		//}
		toc_op(t_list_sssp, 4, tmp_R->num, 0, 88, 0);
#else
		tic_sr(t_list_sssp, 3);
		request_arr_single(gh, S, R);
		toc_op(t_list_sssp, 3, S->num, R->num, 88, 64);

		//q_sort_arr(R, 0, R->num-1);
		//reorder_arr(R, gs);

		tic_sr(t_list_sssp, 4);
		relax_arr_single(gs, predecesor, B, R, delta, &b_total, vidx);
		toc_op(t_list_sssp, 4, R->num, 0, 88, 0);
#endif
		empty_arr(S);
		print_B_map(log, B, iter);
		//printf("iter %d +++++++++++++++++++\n", iter);
		//if(iter==5)
		//	exit(1);
		iter++;
	}
	print_sssp_results(t_list_sssp);
	//DPRINTF(1, "iterations taken %lld\n", iter);	
	//for(i=0;i<gs->v_size;i++)
	//	DPRINTF(1, "vertex %lld %lf\n", (i+1), gs->vet_info[i].weight);
	//DPRINTF(1, "\n");
}

//lint
//check_total(lint *b_total){
//	int i;
//	lint total =0;
//	for(i=0;i<num_threads;i++)
//		total+=b_total[i];
//	return total; 
//}


void recover_arr(p_arr *B, p_vidx vidx, 
		lint iter){
	//this is the recover step
	lint i, j=0;
	//could be done in parallel
	for(i=0; i<B[iter]->idx; i++){
		if(B[iter]->arr[i].int_val != -1){
			tmp_B->arr[j].int_val = B[iter]->arr[i].int_val;
			tmp_B->arr[j].int_val1 = B[iter]->arr[i].int_val1;
			tmp_B->arr[j].double_val = B[iter]->arr[i].double_val;
			lint pos = B[iter]->arr[i].int_val;
			vidx[pos].buck_pos = j;
			j++;	
		}
	}
	B[iter]->idx = j;
	B[iter]->num = j;
	//copy back parallel
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads) shared(B)
#endif
	for(i=0; i<B[iter]->idx; i++){
		B[iter]->arr[i].int_val = tmp_B->arr[i].int_val;
		B[iter]->arr[i].int_val1 = tmp_B->arr[i].int_val1;
		B[iter]->arr[i].double_val = tmp_B->arr[i].double_val;
	}
}

void
relax_arr_single(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx){
	double delta_recip = (double)1.0/delta; 
	lint i, j;
	for(i=0; i<R->idx; i++){
		lint source_idx = R->arr[i].int_val;
		double source_val = gs->vet_info[source_idx].weight;
		lint target_idx = R->arr[i].int_val1;
		//lint edge_idx = find_edge_idx_arr(gs, source_idx, target_idx);
		double target_ori_val = gs->vet_info[target_idx].weight;
		double target_val = R->arr[i].double_val;
		double target_new_val = target_val+source_val;
		//if(target_idx == 22 && iter==5)
		//	printf("source %lld target %lld o_val %lf new_val %lf edge_val %lf source_val %lf\n", source_idx, target_idx, target_ori_val, target_new_val, target_val, source_val);
		if(target_new_val <= target_ori_val){
			gs->vet_info[target_idx].weight = target_new_val;
			lint new_buck;
			if(target_ori_val!=INFINITY){
				lint old_buck = vidx[target_idx].buck_id;
				new_buck = (lint)(target_new_val*delta_recip);	
				if(old_buck!=new_buck){
					if(vidx[target_idx].buck_pos != -1){
						remove_find_arr(B[old_buck], vidx[target_idx].buck_pos);
						B[old_buck]->num--;
					}
					add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, B[new_buck]->idx);
					B[new_buck]->idx++;
					B[new_buck]->num++;
				}
			}
			else{
				//printf("source %lld target %lld o_val %lf new_val %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
				new_buck = (lint)(target_new_val*delta_recip);
				add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, B[new_buck]->idx);
				*(b_total)+=1;
				B[new_buck]->idx++;
				B[new_buck]->num++;
			}
			//change vidx
			vidx[target_idx].buck_id = new_buck;
			vidx[target_idx].buck_pos = B[new_buck]->idx-1;
		}
	}
}

//void
//relax_arr_single(t_csr *gs, 
//		p_list *predecesor, p_arr *B, 
//		p_arr R, double delta, 
//		lint *b_total, p_vidx ;vidx){
//	lint i, j;
//	for(i=0; i<R->idx; i++){
//		lint source_idx = R->arr[i].int_val;
//		double source_val = gs->vet_info[source_idx].weight;
//		lint target_idx = R->arr[i].int_val1;
//		double target_ori_val = gs->vet_info[target_idx].weight;
//		double target_new_val = R->arr[i].double_val;
//		lint pre_target_idx = i==0?-1:R->arr[i-1].int_val1;
//		if(target_new_val <= target_ori_val){
//			gs->vet_info[target_idx].weight = target_new_val;
//			lint new_buck;
//			if(target_ori_val!=INFINITY){
//				lint old_buck = vidx[target_idx].buck_id;
//				new_buck = (lint)(target_new_val/delta);
//				if(old_buck!=new_buck){
//					remove_find_arr(B[old_buck], vidx[target_idx].buck_pos);
//					B[old_buck]->num--;
//					add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, B[new_buck]->idx);
//					B[new_buck]->idx++;
//					B[new_buck]->num++;
//				}
//			}	
//			else{
//				//printf("source %lld target %lld o_val %lf new_val %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
//				new_buck = (lint)(target_new_val/delta);
//				add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, B[new_buck]->idx);
//				*(b_total)+=1;
//				B[new_buck]->idx++;
//				B[new_buck]->num++;
//			}
//			//change vidx
//			vidx[target_idx].buck_id = new_buck;
//			vidx[target_idx].buck_pos = B[new_buck]->idx-1;
//		}
//	}
//}

void			
relax_arr_by_partition(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx){
	int i, j;
	double delta_recip = (double)1.0/delta;
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0;i<num_threads;i++){
		for(j=0;j<buck_size;j++){
			b_size[i][j]=0;
			b_start[i][j]=0;
		}
	}
	//scan to get b_start
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads) shared(r_start, predecesor, gs, B, R, b_total, vidx) private(i, j)
#endif
	for(i=0; i<HASH_BINS; i++){
		sint tid = omp_get_thread_num();
		for(j=bin_idx[i];j<bin_idx[i+1];j++){
			lint source_idx = R->arr[j].int_val;
			double source_val = gs->vet_info[source_idx].weight;
			lint target_idx = R->arr[j].int_val1;
			double target_ori_val = gs->vet_info[target_idx].weight;
			double target_val = R->arr[j].double_val;
			double target_new_val = target_val+source_val;
			R->arr[j].double_val = target_new_val;
			if(target_new_val <= target_ori_val){
				lint new_buck;
				new_buck = (lint)(target_new_val*delta_recip);
				b_size[tid][new_buck]++;
			}
		}
	}
	//compute b_id
	lint b_virtual_idx[num_threads];
	for(j=0;j<buck_size;j++){
		lint sum = B[j]->idx;
		b_start[0][j] = B[j]->idx;
		for(i=1; i<num_threads; i++){
			sum += b_size[i-1][j];
			b_start[i][j] = sum;
		}
		sum += b_size[num_threads-1][j];
		if(sum > g_v_num){
			recover_arr(B, vidx, j);
			lint sum = B[j]->idx;
			b_start[0][j] = B[j]->idx;
			for(i=1; i<num_threads; i++){
				sum += b_size[i-1][j];
				b_start[i][j] = sum;
			}
			sum += b_size[num_threads-1][j];
			B[j]->idx = sum;
		}
		else
			B[j]->idx = sum;
	}
	//this could be parallelized if neccessary
	lint local_b_total[100];
	for(i=0;i<num_threads;i++){
		b_virtual_idx[i] = B[i]->idx;
		local_b_total[i] = 0;
		for(j=0;j<buck_size;j++){
			local_rmv_cnt[i][j] = 0;
			local_add_cnt[i][j] = 0;
			b_virtual_idx[i] += b_size[i][j];
		}
	}
	//major computation
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads) shared(r_start, predecesor, gs, B, R, b_total, vidx) private(i, j)
#endif
	for(i=0; i<HASH_BINS; i++){
		sint tid = omp_get_thread_num();
		for(j=bin_idx[i]; j < bin_idx[i+1]; j++){
			//lint source_idx = R->arr[j].int_val;
			//double source_val = gs->vet_info[source_idx].weight;
			lint target_idx = R->arr[j].int_val1;
			double target_ori_val = gs->vet_info[target_idx].weight;
			//double target_val = R->arr[j].double_val;
			double target_new_val = R->arr[j].double_val;
			//if(is_heavy == TRUE && iter == 6)
			//	printf("source %lld target %lld target_old %lf target_new %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
			if(target_new_val <= target_ori_val){
				if(tid==0)
					tic_sr(t_list_sssp, 5);
				gs->vet_info[target_idx].weight = target_new_val;
				lint new_buck;
				if(target_ori_val!=INFINITY){
					lint old_buck = vidx[target_idx].buck_id;
					new_buck = (lint)(target_new_val*delta_recip);	
					if(vidx[target_idx].buck_pos != -1){
						remove_find_arr(B[old_buck], vidx[target_idx].buck_pos);
						local_rmv_cnt[tid][old_buck]++;
					}
					local_add_cnt[tid][new_buck]++;
					add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, b_start[tid][new_buck]);
				}
				else{
					new_buck = (lint)(target_new_val*delta_recip);
					add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, b_start[tid][new_buck]);
					local_b_total[tid]+=1;
					local_add_cnt[tid][new_buck]++;
				}
				//change vidx
				vidx[target_idx].buck_id = new_buck;
				vidx[target_idx].buck_pos = b_start[tid][new_buck];
				b_start[tid][new_buck]++;
				if(tid==0)
					toc_sr(t_list_sssp, 5);
			}
		}	
	}
	//posterior step
	for(i=0;i<num_threads;i++){
		*(b_total)+=local_b_total[i];
	}
	for(i=0;i<num_threads;i++){
		for(j=0;j<buck_size;j++){
			B[j]->num -= local_rmv_cnt[i][j];
			B[j]->num += local_add_cnt[i][j];
		}
	}
	//for(i=0;i<buck_size;i++)
	//	if(B[i]->num != 0)
	//		printf("iter %lld B[%lld] idx is %lld num is %lld\n", iter, i, B[i]->idx, B[i]->num);
}

void
relax_arr(t_csr *gs, 
		p_list *predecesor, p_arr *B, 
		p_arr R, double delta, 
		lint *b_total, p_vidx vidx){
	int i, j;
	double delta_recip = (double)1.0/delta;
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0;i<num_threads;i++){
		for(j=0;j<buck_size;j++){
			b_size[i][j]=0;
			b_start[i][j]=0;
		}
	}
	//scan to get b_start
#ifdef USE_OMP
#pragma omp parallel num_threads(num_threads) shared(r_start, predecesor, gs, B, R, b_total, vidx) private(i, j)
	{
		sint tid = omp_get_thread_num();
		lint start = tid==0?0:r_start[tid*64];
		lint end = tid==num_threads-1?R->num:r_start[(tid+1)*64];
		//		printf("this is threads %lld start is %lld, end is %lld\n", tid, start, end);
#else
		lint start = 0;
		lint end = 1;
		sint tid=0;
#endif
		for(i=start; i<end; i++){
			lint source_idx = R->arr[i].int_val;
			double source_val = gs->vet_info[source_idx].weight;
			lint target_idx = R->arr[i].int_val1;
			//lint edge_idx = find_edge_idx_arr(gs, source_idx, target_idx);
			//prefetch
			//prefetch_dis=1;
			//int target_next = R->arr[i+prefetch_dis].int_val1;
			//_mm_prefetch((char *)(gs->vet_info+target_next), _MM_HINT_T0);
			double target_ori_val = gs->vet_info[target_idx].weight;
			double target_val = R->arr[i].double_val;
			double target_new_val = target_val+source_val;
			R->arr[i].double_val = target_new_val;
			//double target_new_val = R->arr[i].double_val;
			lint pre_target_idx = i==0?-1:R->arr[i-1].int_val1;
			//printf("new val %lf old val %lf\n", target_new_val, target_ori_val);
			if(target_new_val <= target_ori_val && pre_target_idx != target_idx){
				lint new_buck;
				new_buck = (lint)(target_new_val*delta_recip);
				b_size[tid][new_buck]++;
				//printf("%lld %lld\n", new_buck, tid);
				//if((float)target_new_val==(float)0.3)
				//	printf("source %lld target %lld tid %d value %lf\n", source_idx, target_idx, tid, target_new_val);
			}
			//else
			//	printf("source %lld target %lld old val %lf new val %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
		}
#ifdef USE_OMP
	}
#endif
	//compute b_id
	for(j=0;j<buck_size;j++){
		lint sum = B[j]->idx;
		b_start[0][j] = B[j]->idx;
		for(i=1; i<num_threads; i++){
			sum += b_size[i-1][j];
			b_start[i][j] = sum;
		}
		sum += b_size[num_threads-1][j];
		B[j]->idx = sum;
		if(B[j]->idx > g_v_num){
			//printf("recover here %lld %lld\n", B[j]->idx, j);
			recover_arr(B, vidx, j);
		}
	}
	//if(b_size[0][3]!=0)
	//	for(i=0;i<num_threads;i++)
	//		printf("b start is %lld\n", b_start[i][3]);

	//this could be parallelized if neccessary
	lint local_b_total[100];
	for(i=0;i<num_threads;i++){
		local_b_total[i]=0;
		for(j=0;j<buck_size;j++){
			local_rmv_cnt[i][j]=0;
			B[j]->num+=b_size[i][j];
		}
	}

	//lint m,n;
	//if(iter==1 && inner_i==4){
	//	for(m=iter;m<iter+5;m++){
	//		for(n=0;n<num_threads;n++)
	//			printf("%lld ", b_start[n][m]);
	//		printf("\n");
	//	}
	//}

	//major computation
#ifdef USE_OMP
#pragma omp parallel num_threads(num_threads) shared(r_start, predecesor, gs, B, R, b_total, vidx) private(i, j)
	{
		sint tid = omp_get_thread_num();
		lint start = tid==0?0:r_start[tid*64];
		lint end = tid==num_threads-1?R->num:r_start[(tid+1)*64];
#endif
		for(i=start; i<end; i++){
			lint target_idx = R->arr[i].int_val1;
			double target_ori_val = gs->vet_info[target_idx].weight;
			//double target_val = R->arr[i].double_val;
			double target_new_val = R->arr[i].double_val;
			lint pre_target_idx = i==0?-1:R->arr[i-1].int_val1;
			//if(is_heavy == TRUE && iter == 6)
			//	printf("source %lld target %lld target_old %lf target_new %lf\n", source_idx, target_idx, target_ori_val, target_new_val);
#ifdef REV_STATS
			lint source_idx = R->arr[i].int_val;
			double source_val = gs->vet_info[source_idx].weight;
			lint source_buck = (lint)(source_val*delta_recip);
			lint target_buck = (lint)(target_new_val*delta_recip);
			if(target_new_val>target_ori_val || target_buck==source_buck){
				lint target_ori_buck = (lint)(target_ori_val*delta_recip);
				if(pre_target_idx == target_idx)
					rev_stats[tid][iter][CFLT]++;
				else if(target_ori_buck == target_buck)
					rev_stats[tid][iter][SAME]++;
				else if(target_ori_buck < target_buck)
					rev_stats[tid][iter][BACK]++;
				rev_stats[tid][iter][TTL]++;
			}
#endif

			if(target_idx == 41964672)
				printf("iter %lld inner_iter %d\n", iter, inner_i);
			if(target_new_val <= target_ori_val){
				rev_stats[tid][iter][TTL]++;
				rev_stats[tid][iter][SUC]++;
				if(tid==0)
					tic_sr(t_list_sssp, 5);
				if(pre_target_idx == target_idx)
				{
				//	if(target_new_val==target_ori_val){
				//		//printf("the target vertex is %lld \n", pre_target_idx);
				//		lint bc_id = vidx[target_idx].buck_id;
				//		lint bc_pos = vidx[target_idx].buck_pos;
				//		B[bc_id]->arr[bc_pos].int_val = target_idx;
				//		B[bc_id]->arr[bc_pos].double_val = target_new_val;
				//	}
				}
				else{
					gs->vet_info[target_idx].weight = target_new_val;
					lint new_buck;
					lint old_buck=-1;
					if(target_ori_val!=INFINITY){
						old_buck = vidx[target_idx].buck_id;
						new_buck = (lint)(target_new_val*delta_recip);	
						if(old_buck!=new_buck){
							if(vidx[target_idx].buck_pos != -1){
								remove_find_arr(B[old_buck], vidx[target_idx].buck_pos);
								local_rmv_cnt[tid][old_buck]++;
							}
							add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, b_start[tid][new_buck]);
							//if(target_idx == 145454)
							//	printf("init value\n");
							//printf("%lld \n", b_start[tid][new_buck]);
						}
					}
					else{
						new_buck = (lint)(target_new_val*delta_recip);
						add_arr_two_pos(B[new_buck], target_idx, gs->vet_info[target_idx].weight, b_start[tid][new_buck]);
						local_b_total[tid]+=1;
						//printf("%lld \n", b_start[tid][new_buck]);
					}
					//change vidx
					if(old_buck != new_buck){
						vidx[target_idx].buck_id = new_buck;
						vidx[target_idx].buck_pos = b_start[tid][new_buck];
						b_start[tid][new_buck]++;
					}
				}
				if(tid==0)
					toc_sr(t_list_sssp, 5);
			}
		}	
#ifdef USE_OMP
	}
#endif
	for(i=0;i<num_threads;i++){
		*(b_total)+=local_b_total[i];
	}
	//printf("\n");
	//printf("b total is %lld\n", *b_total);
	//this could be parallelized if neccessary
	for(i=0;i<num_threads;i++){
		for(j=0;j<buck_size;j++){
			B[j]->num -= local_rmv_cnt[i][j];
		}
	}
	//for(i=0;i<buck_size;i++)
	//	if(B[i]->num != 0)
	//		printf("iter %lld B[%lld] idx is %lld num is %lld\n", iter, i, B[i]->idx, B[i]->num);
	//	for(i=0;i<B[2]->idx; i++)
	//		printf("%lld|%lld|%lf\n", B[2]->arr[i].int_val, B[2]->arr[i].int_val1, B[2]->arr[i].double_val);
}


lint
find_edge_idx_arr(t_csr *gs, lint source, 
		lint target){
	lint start = source==0?0:gs->vet_idx[source-1];
	lint end = gs->vet_idx[source];
	lint i;
	for(i=start; i<end; i++){
		if(gs->edge_idx[i]==target)
			return i;
	}
	return -1;
}

//void 
//compute_request_idx(sint thread_id, t_csr *gs, t_arr *arr, FILE* log){
//	lint i, j;
//	r_size[thread_id]=0;
//	for(j=0; j<arr->idx; j++){
//		lint query = arr->arr[j].int_val;
//		lint start = query==0?0:gs->vet_idx[query-1];
//		lint end = gs->vet_idx[query];
//		r_size[thread_id]+=(end-start);
//	}
//	printf("thread %d come here\n", thread_id);
//}

//void 
//reduce_request_idx(){
//	lint sum =0;
//	sint i;
//	for(i=0;i<num_threads; i++){
//		r_start[i]=0;
//	}
//	for(i=0;i<num_threads; i++){
//		sum+=r_size[i];
//		r_size[i]=sum;
//	}
//	for(i=num_threads-1; i>=1; i--){
//		r_start[i]=r_size[i-1];
//	}
//	r_start[0]=0;
//}

void
request_arr(t_csr *gs, p_arr arr, p_arr R){
	lint i, j;
	//scan through to decide b start
	for(i=0;i<num_threads;i++){
		r_size[i*64]=0;
		r_start[i*64]=0;
	}
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0; i<arr->idx; i++){
		sint tid= omp_get_thread_num();
		lint query = arr->arr[i].int_val;
		if(query !=-1){ // this position might be deleted
			//prefetch here
			//prefetch_dis=1;
			//int target_next = gs->vet_idx[query-1+prefetch_dis];
			//_mm_prefetch((char *)(gs->vet_idx + target_next), _MM_HINT_T0);
			lint start = query==0?0:gs->vet_idx[query-1];	
			lint end = gs->vet_idx[query];
			r_size[tid*64]+=(end-start);
		}
	}
	lint sum=0;
	r_start[0]=0;
	for(i=1;i<num_threads;i++){
		sum += r_size[(i-1)*64];
		r_start[i*64]=sum;
	}
	sum+=r_size[(num_threads-1)*64];
	//for(i=0;i<num_threads;i++)
	//	printf("%lld ", r_start[i]);
	//printf("\n");
	//for(i=0;i<num_threads;i++)
	//	printf("%lld|%lld ", r_start[i], gs->edge_idx[i]);
	//printf("\n");
	//start query
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads) private(i)
#endif
	for(j=0; j<arr->idx; j++){
		sint tid= omp_get_thread_num();
		lint query = arr->arr[j].int_val;
		if(query != -1){
			lint start = query==0?0:gs->vet_idx[query-1];	
			lint end = gs->vet_idx[query];
			//prefetch_dis=1;
			//int target_next = gs->vet_idx[query-1+prefetch_dis];
			//_mm_prefetch((char *)(gs->vet_idx + target_next), _MM_HINT_T0);
			//if(query==743787)
			//	printf("heiheiheiheihei\n");
			for(i=start; i<end; i++){
				lint idx = r_start[tid*64];
				lint target = gs->edge_idx[i];
				//double target_val = gs->edge_info[i].edge_weight;
				double target_val = gs->edge_info[i].edge_weight;
				add_arr_three_pos(R, query, target, target_val, idx);
				r_start[tid*64]++;
			}
		}
	}
	R->num = sum;
	R->idx = sum;
	//for(i=0;i<R->idx;i++)
	//	printf("%lld|%lld|%1.1f ", R->arr[i].int_val, R->arr[i].int_val1, R->arr[i].double_val);
	//printf("\n");
}

void
request_arr_single(t_csr *gs, p_arr arr, p_arr R){
	lint i, j;
	R->idx=0;
	R->num=0;
	for(j=0; j<arr->idx; j++){
		lint query = arr->arr[j].int_val;
		if(query != -1){
			lint start = query==0?0:gs->vet_idx[query-1];	
			lint end = gs->vet_idx[query];
			for(i=start; i<end; i++){
				lint idx = R->idx;
				lint target = gs->edge_idx[i];
				double target_val = gs->edge_info[i].edge_weight;
				add_arr_three_pos(R, query, target, target_val, idx);
				R->idx++;
				R->num++;
			}
		}
	}
}

	void 
compute_request_partition(p_arr R)
{
	lint i, j, avg;
	avg = R->num/num_threads + 1;
	for(i=0;i<num_threads;i++)
		r_start[i*64] = 0;
	for(i=1;i<num_threads;i++){
		r_start[i*64] = i*avg>=R->num?R->num:i*avg;
		j= r_start[i*64]-1;
		while(1){
			if(R->arr[j].int_val1 == R->arr[j-1].int_val1){
				j--;
				r_start[i*64]--;
			}
			else
				break;
		}
	}
}

//void
//compute_RR_partition(p_arr RR){
//	int i;
//	lint avg = RR->num/num_threads+1;
//	for(i=0;i<num_threads;i++){
//		rr_start[i]=i*avg;
//		rr_end[i]=(i==num_threads-1?RR->num:(i+1)*avg);
//		if(i>0){
//			lint k=i*avg-1;
//			while(1){
//				if(RR->arr[k].int_val1==RR->arr[k+1].int_val1){
//					k--;
//					rr_start[i]--;
//					rr_end[i-1]--;
//				}
//				else
//					break;
//			}
//		}
//	}
//}

//void
//balance_request(p_arr RR, p_arr R, 
//		sint thread_id){
//	lint i,j;	
//	lint start = rr_start[thread_id];
//	lint end = rr_end[thread_id];
//	for(i=start, j=0; i<end; i++, j++){
//		cpy_selected_arr(RR, &(R->arr[j]), i);	
//		R->idx++;
//		R->num++;
//	}
//}

void
remember_arr_single(p_arr B, p_arr S, 
		lint *b_total, p_vidx vidx){
	lint i,j=S->idx;
	for(i=0;i<B->idx;i++){
		if(B->arr[i].int_val != -1){
			S->arr[j].int_val = B->arr[i].int_val;
			S->arr[j].int_val1 = B->arr[i].int_val1;
			S->arr[j].double_val = B->arr[i].double_val;
			B->arr[i].int_val = -1;
			B->arr[i].int_val1 = -1;
			B->arr[i].double_val = -1;
			lint vidx_idx = S->arr[j].int_val;
			vidx[vidx_idx].buck_id = -1;	
			vidx[vidx_idx].buck_pos = -1;
			j++;
		}
	}
	S->idx += B->num;
	S->num += B->num;	
	(*b_total) -= B->num;
	B->idx = 0;
	B->num = 0;
}

void
remember_arr(p_arr B, p_arr S, 
		lint *b_total, p_vidx vidx){
	lint i, j=0;
	for(i=0;i<num_threads;i++){
		s_start[i*64]=0;
		s_size[i*64]=0;
	}
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0;i<B->idx;i++){
		sint tid = omp_get_thread_num();
		if(B->arr[i].int_val != -1)
			s_size[tid*64]++;
	}

	lint sum = S->idx;
	s_start[0]=S->idx;
	for(i=1;i<num_threads;i++){
		sum+=s_size[(i-1)*64];
		s_start[i*64] = sum;
	}
	//for(i=0;i<num_threads;i++)
	//	printf("%lld ", s_start[i]);
	//printf("\n");

#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0;i<B->idx;i++){
		sint tid = omp_get_thread_num();
		lint sidx = s_start[tid*64];
		if(B->arr[i].int_val != -1){ //might be deleted
			S->arr[sidx].int_val = B->arr[i].int_val;
			S->arr[sidx].int_val1 = B->arr[i].int_val1;
			S->arr[sidx].double_val = B->arr[i].double_val;
			B->arr[i].int_val = -1;
			B->arr[i].int_val1 = -1;
			B->arr[i].double_val = -1;
			lint vidx_idx = S->arr[sidx].int_val;
			vidx[vidx_idx].buck_id = -1;	
			vidx[vidx_idx].buck_pos = -1;
			s_start[tid*64]++;
		}	
	}
	(*b_total)-=B->num;
	S->num = S->num + B->num;
	S->idx = S->idx + B->num;
	//printf("ok here the s size is %lld\n", S->idx);
	//for(i=0;i<S->idx;i++)
	//	printf("%lld ", S->arr[i].int_val);
	B->num = 0;
	B->idx = 0;
}

int
check_B_total_on_iter(p_arr **B, lint iter){
	sint i;
	for(i=0;i<num_threads;i++){
		if(B[i][iter]->num>0)
			return TRUE;	
	}
	return FALSE;
}

	void
set_vidx(p_arr B, p_vidx vidx, sint iter)
{
	lint i;
	for(i=0;i<B->num;i++){
		lint target = B->arr[i].int_val1;	
		vidx[target].buck_id = iter;
		vidx[target].buck_pos = i;
	}
}

void print_B_map(FILE *log, p_arr *B, lint iter){
	lint i, j, k;
	fprintf(log, "============After %lld, the B Map is:================\n", iter);
	for(i=0;i<buck_size; i++){
		if(B[i]->num > 0)
			fprintf(log, "Map on iter %lld\n", i);
		else
			continue;
		for(j=0;j<B[i]->idx;j++){
			fprintf(log, "%lld ", B[i]->arr[j].int_val);
		}
		fprintf(log, "\n");
		fflush(log);
	}
	fflush(log);
}
void
print_sssp_results(p_time_list pl){
	printf("======time for light edge request===========\n");
	printf("cycles %lld, time(seconds) %lf operations %lld bytes %lld BW %lf\n", pl->list[0].cycle, pl->list[0].time, pl->list[0].operations, pl->list[0].bytes, (double)pl->list[0].bytes/pl->list[0].time/1.0e8);
	printf("======time for remember===========\n");
	printf("cycles %lld, time(seconds) %lf operations %lld bytes %lld BW %lf\n", pl->list[1].cycle, pl->list[1].time, pl->list[1].operations, pl->list[1].bytes, (double)pl->list[1].bytes/pl->list[1].time/1.0e8);
	printf("======time for light egde relax===========\n");
	printf("cycles %lld, time(seconds) %lf operations %lld bytes %lld BW %lf\n", pl->list[2].cycle, pl->list[2].time, pl->list[2].operations, pl->list[2].bytes, (double)pl->list[2].bytes/pl->list[2].time/1.0e8);
	printf("======time for heavy edge request===========\n");
	printf("cycles %lld, time(seconds) %lf operations %lld bytes %lld BW %lf\n", pl->list[3].cycle, pl->list[3].time, pl->list[3].operations, pl->list[3].bytes, (double)pl->list[3].bytes/pl->list[3].time/1.0e8);
	printf("======time for heavy egde relax===========\n");
	printf("cycles %lld, time(seconds) %lf operations %lld bytes %lld BW %lf\n", pl->list[4].cycle, pl->list[4].time, pl->list[4].operations, pl->list[4].bytes, (double)pl->list[4].bytes/(pl->list[4].time*1.0e8));
	printf("======time for relax one edge===========\n");
	printf("cycles %lld, time(seconds) %lf\n",  pl->list[5].cycle, pl->list[5].time);
	printf("======time for sorting===========\n");
	printf("cycles %lld, time(seconds) %lf\n",  pl->list[6].cycle, pl->list[6].time);
	lint i,j;
	for(i=0;i<buck_size;i++){
		lint sum_suc=0;
		lint sum_cflt=0;
		lint sum_back=0;
		lint sum_same=0;
		lint sum_ttl=0;
		for(j=0;j<num_threads;j++){
			sum_suc += rev_stats[j][i][SUC];
			sum_cflt += rev_stats[j][i][CFLT];
			sum_back += rev_stats[j][i][BACK];
			sum_same += rev_stats[j][i][SAME];
			sum_ttl += rev_stats[j][i][TTL];
		}
		if(sum_ttl!=0)
			printf("iter %lld: suc: %lld cflt %lld back %lld same %lld ttl %lld\n", i, sum_suc, sum_cflt, sum_back, sum_same, sum_ttl);
	}
}

void 
reorder_arr(p_arr arr, t_csr *gs, double *tmp){
	lint i, j;
#ifdef USE_OMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for(i=0; i<arr->idx; i++){
		//lint source_idx = arr->arr[i].int_val;
		//double source_val = gs->vet_info[source_idx].weight;
		lint target_idx = arr->arr[i].int_val1;
		double target_ori_val = gs->vet_info[target_idx].weight;
		//double target_val = arr->arr[i].double_val;
		double target_new_val = arr->arr[i].double_val;
		//printf("%lf %lf\n", target_new_val, target_ori_val);
		//if(target_new_val <= target_ori_val){
		tmp[i] = target_new_val;
		//}
		//else
		//	tmp[i] = target_ori_val;
	}
	//set partition idx
	lint p_start[num_threads*64];
	lint avg = arr->idx/num_threads +1;
	lint sum =0;
	p_start[0] = 0;
	for(i=1;i<num_threads;i++){
		sum += avg;
		p_start[i*64] = sum;
	}
	for(i=1;i<num_threads;i++){
		while(1){
			lint pre = arr->arr[p_start[i*64]-1].int_val1;
			lint now = arr->arr[p_start[i*64]].int_val1;
			if(pre ==  now)
				p_start[i*64]--;
			else
				break;
		}
	}
	//for(i=0;i<num_threads;i++)
	//	printf("%lld ", p_start[i*64]);
	//printf("\n");
#ifdef USE_OMP
#pragma omp parallel num_threads(num_threads) private(i)
	{
#endif
		sint tid = omp_get_thread_num();
		lint end = tid==(num_threads-1)?arr->idx-1:p_start[(tid+1)*64]-1;
		lint start = p_start[tid*64]+1;
		for(i=end; i>=start; i--){
			if(arr->arr[i-1].int_val1 == arr->arr[i].int_val1 && tmp[i-1] > tmp[i]){
				lint tmp_int = arr->arr[i].int_val;
				lint tmp_int1 = arr->arr[i].int_val1;
				double tmp_double = arr->arr[i].double_val;
				arr->arr[i].int_val = arr->arr[i-1].int_val;
				arr->arr[i].int_val1 = arr->arr[i-1].int_val1;
				arr->arr[i].double_val = arr->arr[i-1].double_val;
				arr->arr[i-1].int_val = tmp_int;
				arr->arr[i-1].int_val1 = tmp_int1;
				arr->arr[i-1].double_val = tmp_double;
				tmp_double = tmp[i];
				tmp[i] = tmp[i-1];
				tmp[i-1] = tmp_double;
			}
		}
#ifdef USE_OMP
	}
#endif
	//if(iter==1 && inner_i==4)
	//for(i=0;i<arr->idx;i++){
	//	printf("%lld|%lf ", arr->arr[i].int_val1, tmp[i]);
	//}
	//printf("\n");
}
