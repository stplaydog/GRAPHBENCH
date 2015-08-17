#include "graph.h"
#include <mpi.h>

static clockid_t clockid;

struct time_elem
{
	unsigned long long start;
	unsigned long long end;
	unsigned long long cycle;
	unsigned long long operations;
	unsigned long long bytes;
	double time;
};
typedef struct time_elem t_time_elem;
typedef struct time_elem *p_time_elem;

struct time_list
{
	t_time_elem *list;	
	int num;
};
typedef struct time_list t_time_list;
typedef struct time_list *p_time_list;

void toc_ser (p_time_list pl, int id);
void tic_ser (p_time_list pl, int id);
double timer_ser ();

void toc_mpi (p_time_list pl, int id);
void tic_mpi (p_time_list pl, int id);
double timer_mpi ();

void toc_op(p_time_list pl, int id, lint operations, lint operations1, lint byte_per_op, lint byte_per_op1);

