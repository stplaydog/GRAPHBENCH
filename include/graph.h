#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <immintrin.h>
#include "barrier.h"

#ifndef _H_GRAPH
#define _H_GRAPH

#define FALSE 0
#define TRUE 1
#define RDM_JMP 0.3
#define THRESH 0.000001
#define HASH_BASE 16769023
#define TYPE_INT 0
#define TYPE_DOUBLE 1
#define RECV 0
#define SEND 1
#define RECV_FIRST 0
#define SEND_FIRST 1
#define MAX_MSG_SZ 100000
#define MASTER 0
#ifdef USE_GRAPH500
#define OFFSET 0
#else
#define OFFSET 1
#endif
#define TYPE_LIGHT 0
#define TYPE_HEAVY 1
#define UNDEFINED -1
#define BUFFSIZE 10000;

#define MPI_BLOCK

#ifdef USE_DEBUG
#define _DEBUG_LEVEL_    3     //print debug info 
#else
#define _DEBUG_LEVEL_    1     //not print debug info
#endif     

#define max(a, b) ((a) > (b) ? (a) : (b))

#define min(a, b) ((a) < (b) ? (a) : (b))

#if ( _DEBUG_LEVEL_ == -1 )
#define DPRINTF( level, fmt, ... )        {}
#else
#define DPRINTF( level, fmt, args... )              \
        do                                                           \
        {                                                             \
            if ( (unsigned)(level) <= _DEBUG_LEVEL_ ) \
            {                                                         \
                fprintf( stdout, fmt, ##args );             \
                fflush( stdout );                                 \
            }                                                         \
        } while ( 0 )
#endif


#define ERROR_EXIT(fmt, args... )                    \
    do                                                           \
    {                                                             \
        fprintf (stderr, "Error: ");                          \
        fprintf( stderr, fmt, ##args );                  \
        fflush( stderr );                                      \
        exit (0);                                                \
    } while ( 0 )

#define ERROR_PRINT() {printf("Error on (%d) line in (%s) file\n", __LINE__, __FILE__); exit(123);}

typedef int sint;
typedef long long int lint;

extern sint num_threads;
extern sint num_nodes;
extern char *par_extend;
extern char *data_extend;
extern char *idx_extend;
extern char *base_dir;
extern char *nodes_extension;
extern lint g_v_num;
extern int bin;
extern int buck_size;
extern int buff_size;
extern lint prefetch_dis;
extern sint use_partition;
//extern t_time_list *t_list;


struct vet_info
{
	double weight; //could be uesed as distance
	double wb;
	double recip;
	lint vet_deg;
};
typedef struct vet_info t_vi;
typedef struct vet_info *p_vi;

struct edge_info
{
	lint node_id;
	double edge_weight; //this is for the single source shortest path
};
typedef struct edge_info t_e;
typedef struct edge_info *p_e;

struct CSR
{
	lint *vet_idx;
	lint *edge_idx;
	t_vi *vet_info;	
	t_e *edge_info;

	lint v_size;
	lint e_size;
	lint offset;
	sint csr_type; //could be send or recv
	lint global_v_size;
};
typedef struct CSR t_csr;
typedef struct CSR *p_csr;

struct MSG
{
	double **double_send;
	double **double_recv;
	lint **int_send;
	lint **int_recv;
	lint *idx_send;
	lint *idx_recv;
};
typedef struct MSG t_msg;
typedef struct MSG *p_msg;

struct msg_buff
{
	sint size;
	sint **buff_idx;
	double ***double_buffer;
	lint ***int_buffer;
};
typedef struct msg_buff t_buff;
typedef struct msg_buff *p_buff;


#define HASH_BINS (256)
#define NUMENTRIESPERTABLE (HASH_BINS)
static TREE_BARRIER_TYPE tree_barrier;
#define __BARRIER__ TREE_BARRIER(&tree_barrier, threadid, nthreads);

#endif
