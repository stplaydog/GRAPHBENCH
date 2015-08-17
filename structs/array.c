#include "array.h"
#include "barrier.h"

lint 
bin_search_fuzzy(p_arr arr, lint val, lint start, lint end);

void
create_arr(p_arr arr, lint size){
	lint i;
	arr->arr = (t_arr_elem*)malloc(sizeof(t_arr_elem)*size);	
	arr->size = size;
	arr->num = 0;
	arr->idx = 0;
	for(i=0; i<size; i++){
		arr->arr[i].int_val =-1;
		arr->arr[i].int_val1 =-1;
		arr->arr[i].double_val =-1;
	}
}

void 
add_arr_one(p_arr arr, lint val1){
	lint idx = arr->idx++;
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = -1;	
	arr->arr[idx].double_val = -1;	
	arr->num++;
}

void 
add_arr_two(p_arr arr, lint val1, double dval1){
	lint idx = arr->idx++;
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = -1;	
	arr->arr[idx].double_val = dval1;	
	arr->num++;
}

void 
add_arr_three(p_arr arr, lint val1, lint val2, double dval1){
	lint idx = arr->idx++;
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = val2;	
	arr->arr[idx].double_val = dval1;	
	arr->num++;
}

void 
add_arr_one_pos(p_arr arr, lint val1, 
		lint idx){
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = -1;	
	arr->arr[idx].double_val = -1;	
}

void 
add_arr_two_pos(p_arr arr, lint val1, 
		double dval1, lint idx){
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = -1;	
	arr->arr[idx].double_val = dval1;	
}

void 
add_arr_three_pos(p_arr arr, lint val1, 
		lint val2, double dval1,
		lint idx){
	arr->arr[idx].int_val = val1;	
	arr->arr[idx].int_val1 = val2;	
	arr->arr[idx].double_val = dval1;	
}

void
pop_elem_arr(p_arr arr, p_arr_elem elem){
	int idx = arr->idx -1;
	(*elem).int_val = arr->arr[idx].int_val;	
	(*elem).int_val1 = arr->arr[idx].int_val1;	
	(*elem).double_val = arr->arr[idx].double_val;	
	arr->idx -= 1;
	arr->num--;
}

void
cpy_selected_arr(p_arr arr, p_arr_elem elem, lint idx)
{
	(*elem).int_val = arr->arr[idx].int_val;	
	(*elem).int_val1 = arr->arr[idx].int_val1;	
	(*elem).double_val = arr->arr[idx].double_val;	
}

void 
remove_find_arr(p_arr arr, lint idx){
	arr->arr[idx].int_val = -1;	
	arr->arr[idx].int_val1 = -1;	
	arr->arr[idx].double_val = -1;
}

void 
empty_arr(p_arr arr){
	int i;
#ifdef USE_OMP
#pragma omp parallel
#endif
	for(i=0; i<arr->idx; i++)
		remove_find_arr(arr, i);
	arr->num=0;
	arr->idx=0;
}

void q_sort_arr(p_arr arr, lint left, lint right){
	lint pivot_int, pivot_int1, l_hold, r_hold;
	double pivot_double;
	l_hold = left;
	r_hold = right;
	pivot_int = arr->arr[left].int_val;
	pivot_int1 = arr->arr[left].int_val1;
	pivot_double = arr->arr[left].double_val;
	while (left < right)
	{
		while ((arr->arr[right].int_val1 >= pivot_int1) && (left < right))
			right--;
		if (left != right)
		{
			arr->arr[left].int_val = arr->arr[right].int_val;
			arr->arr[left].int_val1 = arr->arr[right].int_val1;
			arr->arr[left].double_val = arr->arr[right].double_val;
			left++;
		}
		while ((arr->arr[left].int_val1 <= pivot_int1) && (left < right))
			left++;
		if (left != right)
		{
			arr->arr[right].int_val = arr->arr[left].int_val;
			arr->arr[right].int_val1 = arr->arr[left].int_val1;
			arr->arr[right].double_val = arr->arr[left].double_val;
			right--;
		}
	}
	arr->arr[left].int_val = pivot_int;
	arr->arr[left].int_val1 = pivot_int1;
	arr->arr[left].double_val = pivot_double;
	pivot_int1 = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot_int1)
		q_sort_arr(arr, left, pivot_int1-1);
	if (right > pivot_int1)
		q_sort_arr(arr, pivot_int1+1, right);
}


void
partition(p_arr arr, p_arr tmp, lint * bin_idx)
{
  int N = arr->idx;
  int* histogram = (int*) _mm_malloc(sizeof(int)*NUMENTRIESPERTABLE*num_threads, 64);
  int* position = (int*) _mm_malloc(sizeof(int)*N, 64);
  int right_shift = (int)(log2(g_v_num))-8;
  
  #pragma omp parallel num_threads(num_threads)
  {
    int i, t; 
    int threadid = omp_get_thread_num();
    int p_per_thread = N/num_threads;
    if((p_per_thread * num_threads) != N) p_per_thread++;
    int start_n = p_per_thread*threadid;
    int end_n = start_n + p_per_thread;
    if (end_n > N) end_n = N;

    memset(histogram + threadid*NUMENTRIESPERTABLE, 0, sizeof(int)*NUMENTRIESPERTABLE);

    for (i = start_n; i < end_n; i++) {
      int finalHash = arr->arr[i].int_val1 >> right_shift;
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
      if(current_sum != N) ERROR_PRINT();
    }
//__BARRIER__;
TREE_BARRIER(&tree_barrier, threadid, num_threads);

    for (i = start_n; i < end_n; i++) {
      int finalHash = arr->arr[i].int_val1 >> right_shift;
      int newIndex = histogram[threadid*NUMENTRIESPERTABLE + finalHash] + position[i];
      tmp->arr[newIndex] = arr->arr[i];
    }

  }

  _mm_free(histogram);
  _mm_free(position);
}

void 
merge_sort(p_arr arr, p_arr tmp){
	lint i, j;
	//copy from
	if(tmp->num < 1024){
		arr->idx = tmp->idx;
		arr->num = tmp->num;
		q_sort_arr(tmp, 0, arr->idx -1);
		for(i=0;i<arr->idx;i++){
			arr->arr[i].int_val = tmp->arr[i].int_val;
			arr->arr[i].int_val1 = tmp->arr[i].int_val1;
			arr->arr[i].double_val = tmp->arr[i].double_val;
		}
		return;
	}
	arr->idx = tmp->idx;
	arr->num = tmp->num;
	//#ifdef USE_OMP
	//#pragma omp parallel for num_threads(num_threads)
	//#endif
	//	for(i=0;i<arr->idx;i++){
	//		tmp->arr[i].int_val = arr->arr[i].int_val;
	//		tmp->arr[i].int_val1 = arr->arr[i].int_val1;
	//		tmp->arr[i].double_val = arr->arr[i].double_val;
	//	}
	//partial sort
	lint avg = arr->num/num_threads+1;
	lint tmp_s[num_threads];
#ifdef USE_OMP
#pragma omp parallel num_threads(num_threads) shared(avg, tmp, tmp_s)
	{
#endif
		sint tid = omp_get_thread_num();
		lint start = tid*avg;
		tmp_s[tid] = start;
		lint end = tid==(num_threads-1)?arr->idx-1:(tid+1)*avg-1;
		q_sort_arr(tmp, start, end);
#ifdef USE_OMP
	}
#endif
	//printf("===============================================\n");
	//for(j=0;j<num_threads;j++){
	//	lint start = tmp_s[j];
	//	lint end = j==(num_threads-1)?arr->idx:tmp_s[j+1];
	//	printf("thread %lld\n", j);
	//	for(i=start; i<end ;i++)
	//		printf("%lld ", tmp->arr[i].int_val1);
	//	printf("\n");
	//}
	//exit(1);
#ifdef USE_PMG
	//decide read position
	lint **p_start = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64); 
	lint **p_end = (lint**)_mm_malloc(sizeof(lint*)*num_threads, 64); 
	for(i=0;i<num_threads;i++){
		p_start[i]=(lint*)_mm_malloc(sizeof(lint*)*num_threads, 64);
		p_end[i]=(lint*)_mm_malloc(sizeof(lint*)*num_threads, 64);
	}
	//decide position of each thread
	avg = tmp_s[1]/num_threads+1;
	for(i=0; i<num_threads; i++){
		p_start[i][0] = i*avg;
		p_end[i][0] = (i==num_threads-1)?tmp_s[1]:(i+1)*avg;
		lint end_val = tmp->arr[p_end[i][0]-1].int_val1;
		for(j=1; j<num_threads; j++){
			lint start = tmp_s[j];
			lint end = j==(num_threads-1)?arr->idx-1:tmp_s[j+1]-1;
			lint pos = bin_search_fuzzy(tmp, end_val, start, end);
			while(1 && pos < end+1){
				if(tmp->arr[pos].int_val1>end_val)
					pos--;
				else{
					pos++;
					break;
				}
			}
			while(1 && pos < end+1){
				if(tmp->arr[pos].int_val1<=end_val)
					pos++;
				else
					break;
			}
			p_start[i][j] = i==0?tmp_s[j]:p_end[i-1][j];
			p_end[i][j] = i==(num_threads-1)?(end+1):pos;
		}
	}
	//for(i=0;i<num_threads;i++){
	//	for(j=0;j<num_threads;j++)
	//		printf("%lld|%lld ", p_start[i][j], p_end[i][j]);
	//	printf("\n");
	//}
	//decide write position
	lint arr_start[num_threads];
	lint arr_sum[num_threads];
	for(i=0;i<num_threads;i++){
		arr_sum[i] = 0;
		for(j=0;j<num_threads;j++)
			arr_sum[i] += (p_end[i][j]-p_start[i][j]);
	}
	lint sum=0;
	arr_start[0]=0;
	for(i=1;i<num_threads;i++){
		sum+=arr_sum[i-1];
		arr_start[i]=sum;
	}
	//for(i=0;i<num_threads;i++){
	//	lint start = tmp_s[i];
	//	lint end = i==(num_threads-1)?arr->idx:tmp_s[i+1];
	//	for(j=start;j<end;j++)
	//		printf("%lld:%lld ", j, tmp->arr[j].int_val1);
	//	printf("\n");
	//}
	//printf("\n");
	//for(i=0;i<num_threads;i++)
	//	printf("%lld ", arr_start[i]);
	//printf("\n");

	//parallel merge
#pragma omp parallel num_threads(num_threads) private(i, j) shared(arr_start, p_start, p_end)
	{
		sint tid = omp_get_thread_num();
		lint s[num_threads];
		lint t[num_threads];
		for(i=0;i<num_threads;i++){
			s[i] = p_start[tid][i];
			t[i] = p_end[tid][i];
		}
		//if(tid == 1)
		//	for(i=0;i<num_threads;i++)
		//		printf("%lld|%lld ", s[i], t[i]);
		//printf("\n");
		i=arr_start[tid];
		lint end = tid==(num_threads-1)?arr->idx:arr_start[tid+1];
		//for(j=0;j<num_threads;j++){
		//	lint k;
		//	for(k=s[j];k<t[j];k++){
		//		arr->arr[i].int_val = tmp->arr[k].int_val;
		//		arr->arr[i].int_val1 = tmp->arr[k].int_val1;
		//		arr->arr[i].double_val = tmp->arr[k].double_val;
		//		i++;
		//	}
		//}
		//q_sort_arr(arr, arr_start[tid], end-1);

		while(i< end){
			lint min = 10000000000000;
			lint min_id =0;
			for(j=0;j<num_threads;j++){
				lint id = s[j];
				lint term = t[j];
				if(id>=term)
					continue; //there is no more work for this thread
				if(tmp->arr[id].int_val1 < min){
					min = tmp->arr[id].int_val1;
					min_id = j;
				}
			}
			lint idx = s[min_id];
			arr->arr[i].int_val = tmp->arr[idx].int_val;
			arr->arr[i].int_val1 = tmp->arr[idx].int_val1;
			arr->arr[i].double_val = tmp->arr[idx].double_val;
			s[min_id]+=1;
			i+=1;
		}
	}
	//for(i=1;i<arr->idx;i++){
	//	printf("%lld:%lld ", i, arr->arr[i].int_val1);
	//	if(arr->arr[i].int_val1 < arr->arr[i-1].int_val1)
	//		printf("error ");
	//}
	//printf("\n");
	//exit(1);
#else
	//copy back
	lint s[num_threads];
	for(i=0;i<num_threads;i++)
		s[i] = tmp_s[i];
	i=0;
	while(i< arr->num){
		lint min = 10000000000000;
		lint min_id =0;
		for(j=0;j<num_threads;j++){
			lint id = s[j];
			lint term = j==(num_threads-1)?arr->idx-1:tmp_s[j+1]-1;
			if(id>term)
				continue; //there is no more work for this thread
			if(tmp->arr[id].int_val1 < min){
				min = tmp->arr[id].int_val1;
				min_id = j;
			}
		}
		lint idx = s[min_id];
		arr->arr[i].int_val = tmp->arr[idx].int_val;
		arr->arr[i].int_val1 = tmp->arr[idx].int_val1;
		arr->arr[i].double_val = tmp->arr[idx].double_val;
		//printf("%lld ", arr->arr[i].int_val1);
		s[min_id]++;
		i++;
	}
	//for(i=1; i< arr->num; i++){
	//	printf("%lld ", arr->arr[i].int_val1);
	//	if(arr->arr[i].int_val1 < arr->arr[i-1].int_val1)
	//		printf("error\n");
	//}
	//printf("\n");
	//exit(1);
#endif
}

lint bin_search_fuzzy(p_arr arr, lint val, lint start, lint end){
	lint i;
	lint mid = start + (end -start)/2;
	if(start==end){
		return start;
	}
	else if(start==(end-1)){
		if (val<=arr->arr[start].int_val1)
			return start;
		else
			return end;
	}
	else if(val<=arr->arr[mid].int_val1 && val>arr->arr[mid-1].int_val1)
		return mid;
	else if(val>arr->arr[mid].int_val1)
		return bin_search_fuzzy(arr, val, mid, end);
	else if(val<=arr->arr[mid-1].int_val1)
		return bin_search_fuzzy(arr, val, start, mid);
	//while(1){
	//	if(arr->arr[mid+1].int_val1 <=val)
	//		mid += 1;
	//	else
	//		break;
	//}
	return mid;
}
