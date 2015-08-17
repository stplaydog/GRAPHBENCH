#include "timer.h"

unsigned long long
timer_sr ()
{
  return _rdtsc();
}

void init_timer(p_time_list pl, int size){
	pl->num = size;
	int i;
	for(i=0; i<size; i++){
		pl->list[i].start = 0;
		pl->list[i].end = 0;
		pl->list[i].cycle = 0;
		pl->list[i].time = 0;
		pl->list[i].operations = 0;
		pl->list[i].bytes = 0;
	}	
}


void
tic_sr (p_time_list pl, int id)
{
  pl->list[id].start = timer_sr ();
}

void
toc_sr (p_time_list pl, int id)
{
  pl->list[id].end = timer_sr ();
  pl->list[id].cycle += (pl->list[id].end - pl->list[id].start);
  pl->list[id].time += (double)(pl->list[id].end - pl->list[id].start)/(double)2.7e9; 
}

void toc_op(p_time_list pl, int id, lint operations, lint operations1, lint byte_per_op, lint byte_per_op1){
  pl->list[id].end = timer_sr ();
  pl->list[id].cycle += (pl->list[id].end - pl->list[id].start);
  pl->list[id].time += (double)(pl->list[id].end - pl->list[id].start)/(double)2.7e9; 
  pl->list[id].operations += operations + operations1; 
  pl->list[id].bytes += byte_per_op*operations+ byte_per_op1*operations1; 
}

double
timer_mpi ()
{
  return MPI_Wtime();
}


void
tic_mpi (p_time_list pl, int id)
{
  pl->list[id].start = timer_mpi ();
}

void
toc_mpi (p_time_list pl, int id)
{
  pl->list[id].end = timer_mpi ();
  pl->list[id].time = pl->list[id].end-pl->list[id].start;
}

