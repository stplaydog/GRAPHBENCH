#include "utils.h"
#include <stdarg.h>


lint hash(lint a)
{
	return a%HASH_BASE;
}

void q_sort(lint *numbers, lint left, lint right)
{
	lint pivot, l_hold, r_hold;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right))
			right--;
		if (left != right)
		{
			numbers[left] = numbers[right];
			left++;
		}
		while ((numbers[left] <= pivot) && (left < right))
			left++;
		if (left != right)
		{
			numbers[right] = numbers[left];
			right--;
		}
	}
	numbers[left] = pivot;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort(numbers, left, pivot-1);
	if (right > pivot)
		q_sort(numbers, pivot+1, right);
}

void q_sort_two(lint *key, lint *val, lint left, lint right)
{
	lint pivot, l_hold, r_hold, val_pivot;

	l_hold = left;
	r_hold = right;
	pivot = key[left];
	val_pivot = val[left];
	while (left < right)
	{
		while ((key[right] >= pivot) && (left < right))
			right--;
		if (left != right)
		{
			key[left] = key[right];
			val[left] = val[right];
			left++;
		}
		while ((key[left] <= pivot) && (left < right))
			left++;
		if (left != right)
		{
			key[right] = key[left];
			val[right] = val[left];
			right--;
		}
	}
	key[left] = pivot;
	val[left] = val_pivot;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort_two(key, val, left, pivot-1);
	if (right > pivot)
		q_sort_two(key, val, pivot+1, right);
}

void q_sort_double(double *numbers, lint *idx, lint left, lint right)
{
	lint pivot_idx, l_hold, r_hold;
	double pivot;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	pivot_idx = idx[left];
	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right))
			right--;
		if (left != right)
		{
			numbers[left] = numbers[right];
			idx[left] = idx[right];
			left++;
		}
		while ((numbers[left] <= pivot) && (left < right))
			left++;
		if (left != right)
		{
			numbers[right] = numbers[left];
			idx[right] = idx[left];
			right--;
		}
	}
	numbers[left] = pivot;
	idx[left] = pivot_idx;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort_double(numbers, idx, left, pivot-1);
	if (right > pivot)
		q_sort_double(numbers, idx, pivot+1, right);
}

void q_sort_int(lint *index, double *numbers, lint left, lint right)
{
	lint pivot, l_hold, r_hold;
	double pivot_numbers;

	l_hold = left;
	r_hold = right;
	pivot = index[left];
	pivot_numbers = numbers[left];
	while (left < right)
	{
		while ((index[right] >= pivot) && (left < right))
			right--;
		if (left != right)
		{
			index[left] = index[right];
			numbers[left] = numbers[right];
			left++;
		}
		while ((index[left] <= pivot) && (left < right))
			left++;
		if (left != right)
		{
			index[right] = index[left];
			numbers[right] = numbers[left];
			right--;
		}
	}
	index[left] = pivot;
	numbers[left] = pivot_numbers;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort_int(index, numbers, left, pivot-1);
	if (right > pivot)
		q_sort_int(index, numbers, pivot+1, right);
}

void print_1d_array(void *array, lint len, sint type){
	lint i;
	if(type==TYPE_INT){
		for(i=0;i<len;i++)
			printf("%lld ", ((lint*)array)[i]);
	}
	else if(type==TYPE_DOUBLE){
		for(i=0;i<len;i++)
			printf("%f ", ((double*)array)[i]);
	}
	printf("\n");
}

void print_2d_array(void **array, lint len, lint len2, sint type){
	lint i,j;
	if(type==TYPE_INT){
		for(i=0;i<len;i++){
			for(j=0;j<len2;j++)
				printf("%lld ", ((lint*)array[i])[j]);
			printf("\n");
		}
	}
	else if(type==TYPE_DOUBLE){
		for(i=0;i<len;i++){
			for(j=0;j<len2;j++)
				printf("%f ", ((double*)array[i])[j]);
			printf("\n");
		}
	}
	printf("\n");
}

void max1d(void * arr, lint len, sint type, void *result){
	lint i;
	if(type==TYPE_INT){
		for(i=0;i<len;i++)
			if(((lint*)arr)[i]>*((lint*)result))
				*((lint*)result) = ((lint*)arr)[i];
	}
	else if(type==TYPE_DOUBLE){
		for(i=0;i<len;i++)
			if(((double*)arr)[i]>*((double*)result))
				*((double*)result) = ((double*)arr)[i];
	}
}

void max2d(void **arr, lint len, lint len2, sint type, void *result){
	lint i,j;
	if(type==TYPE_INT){
		for(i=0;i<len;i++)
			for(j=0;j<len2;j++)
				if(((lint*)arr[i])[j]>*((lint*)result))
					*((lint*)result) = ((lint*)arr[i])[j];
	}
	else if(type==TYPE_DOUBLE){
		for(i=0;i<len;i++)
			for(j=0;j<len2;j++)
				if(((double*)arr[i])[j]>*((double*)result))
					*((double*)result) = ((double*)arr[i])[j];
	}
}

lint bin_search(lint *distrib, lint val, lint start, lint end){
	lint mid = start + (end -start)/2;
	if(start==end)
		return start;
	else if(start==(end-1)){
		if (val<=distrib[start])
			return start;
		else
			return end;
	}
	else if(val<=distrib[mid] && val>distrib[mid-1])
		return mid;
	else if(val>distrib[mid])
		return bin_search(distrib, val, mid, end);
	else if(val<=distrib[mid-1])
		return bin_search(distrib, val, start, mid);
	return -1;
}


void concate(char *buf, sint count, ...)
{
   va_list list;
   sint j = 0;
   va_start(list, count); 
   strcpy(buf, va_arg(list, char*));
   for(j=1; j<count; j++)
   {
	   strcat(buf, va_arg(list, char*));
   }
   va_end(list);
}

int str_compare(char *str1, char* str2){
	int i=0;
	while(1){
		if(str1[i]==str2[i])
			;
		else if(str1[i]!=str2[i])
			return FALSE;
		else if(str1[i]=='\0' && str2[i]=='\0')
			break;
		i++;
	}
	return TRUE;
}
