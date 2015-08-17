#pragma once
#include "graph.h"

#ifndef _H_UTILS
#define _H_UTILS 
lint hash(lint a);
lint bin_search_bucket(lint *distrib, lint val, lint start, lint end);
void q_sort(lint *numbers, lint left, lint right);
void q_sort_two(lint *key, lint* val, lint left, lint right);
void q_sort_double(double *numbers, lint *idx, lint left, lint right);
void q_sort_int(lint *index, double *numbers, lint left, lint right);
void print_1d_array(void *array, lint len, sint type);
void print_2d_array(void **array, lint len, lint len2, sint type);
void max1d(void * arr, lint len, sint type, void *result);
void max2d(void **arr, lint len, lint len2, sint type, void *result);
void concate(char *buf, sint count, ...);

#endif
