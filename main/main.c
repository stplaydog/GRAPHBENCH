#include "./graph.h"
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include "timer.h"
#include "list.h"

sint num_threads = 1;
sint num_nodes = 1;
int is_bin=1;
lint g_v_num=1;
int g_e_num=1;
char *par_extend = "parex/";
char *data_extend = "datex/";
char *idx_extend = "idxex/";
char *base_dir = "/home/nsatish/panfs/pagerank/posdata/rmat/";
char *nodes_extension;
int bin=0;
int buck_size=10;
sint use_partition = FALSE;
int buff_size = 10000;
lint prefetch_dis = 10;
//t_time_list *t_list;

int main(int argc, char *argv[])
{
	sint c;
	sint adj = 0;
	sint csr = 0;
	sint chunk =0;
	sint pchunk =0;
	sint bchunk =0;
	sint mpi =0;
	sint delta_stepping = 0;
	sint bc=0;
	sint help=0;
	sint scale;
	sint t2b = 0;
	sint reduce = 0;
	char *file;
	double time=0;
	double start=0, end=0;
	//t_list = (t_time_list*)malloc(sizeof(t_time_list));
	//t_list->list = (t_time_elem*)malloc(sizeof(t_time_elem)*10);

	MPI_Init(&argc, &argv);	

	while ((c = getopt (argc, argv, "tregdyijlakhcp:f:b:n:s:")) != -1)
		switch (c)
		{
			case 'a':
				adj = 1;
				break;
			case 'b':
				base_dir = optarg;
				break;
			case 'c':
				csr = 1;
				break;
			case 'd':
				delta_stepping = 1;
				break;
			case 'e':
				bc = 1;
				break;
			case 'f':
				file = optarg;
				break;
			case 'g':
				pchunk = 1;
				break;
			case 'h':
				help = 1;
				break;
			case 'i':
				mpi = 1;
				break;
			case 'j':
				t2b = 1;
				break;
			case 'l':
				bin = 1;
				break;
			case 'k':
				chunk = 1;
				break;
			case 'n':
				num_nodes = atoi(optarg);
				nodes_extension = (char*)malloc(sizeof(char)*100);
				snprintf(nodes_extension, 10, "%d", num_nodes);
				break;
			case 'p':
				num_threads = atoi(optarg);
				break;
			case 'r':
				reduce = 1;
				break;
			case 's':
				scale = atoi(optarg);
				g_v_num = pow(2, scale);
				g_e_num = g_v_num*16;
				break;
			case 't':
				use_partition = TRUE;
				break;
			case 'y':
				bchunk = 1;
				break;
			case '?':
				fprintf (stderr,
						"Unknown option character `\\x%x'.\n",
						optopt);
				return 1;
			default:
				abort ();
		}
	if(chunk ==1){
		t_csr * gs = (t_csr*)malloc(sizeof(t_csr));
		scan_csr_idx(gs, file, base_dir, SEND);
		int i;
		for(i=0;i<gs->v_size;i++){
			DPRINTF(1, "%d|%lld ", i, gs->vet_info[i].vet_deg);
		}
		DPRINTF(1, "\n");
		//DPRINTF(1, "finished reading recv graph files!\n");
		t_csr * gr = (t_csr*)malloc(sizeof(t_csr));
		scan_csr_idx(gr, file, base_dir, RECV);
		//DPRINTF(1, "finished reading send graph files!\n");
		lint *distrib = (lint*)malloc(sizeof(lint)*num_nodes);
		lint *send_size = (lint*)malloc(sizeof(lint)*num_nodes);
		lint *recv_size = (lint*)malloc(sizeof(lint)*num_nodes);
		comp_distrib_csr(gs, gr, distrib, send_size, recv_size, SEND);
		write_recip(file, base_dir, distrib, gs);
		chunk_files(file, base_dir, distrib, send_size, recv_size);
		return 0;
	}
	if(t2b == 1){
		sint myrank, size;
		MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
		MPI_Comm_size (MPI_COMM_WORLD, &size);
		lint *global_info = (lint*)malloc(sizeof(lint)*(num_nodes+2));
		printf("I am node %d, I am doing chunking now\n", myrank);
		chunk_files_physical_to_bin(file, base_dir, global_info, myrank);
		MPI_Finalize();
	}
	if(pchunk==1){
		printf("here\n");
		//int *ac = (int*)malloc(sizeof(int*));
		//char ***av = (char***)malloc(sizeof(char***));
		//tic_mpi(t_list, 1);
		sint myrank, size;
		MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
		MPI_Comm_size (MPI_COMM_WORLD, &size);
		lint *global_info = (lint*)malloc(sizeof(lint)*(num_nodes+2));
		printf("start physical chunking\n");
		chunk_files_physical(file, base_dir, global_info, myrank);
		printf("end physical chunking\n");
		MPI_Bcast(global_info, (num_nodes+2), MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		lint *distrib = (lint*)malloc(sizeof(lint)*num_nodes);
		lint *send_size = (lint*)malloc(sizeof(lint)*num_nodes);
		lint *recv_size = (lint*)malloc(sizeof(lint)*num_nodes);
		comp_distrib_csr_par(file, base_dir, distrib, send_size, recv_size,global_info,  SEND, myrank);
		chunk_files_par(file, base_dir, distrib, send_size, recv_size, SEND, global_info, myrank);
		chunk_files_par(file, base_dir, distrib, send_size, recv_size, RECV, global_info, myrank);
		MPI_Barrier(MPI_COMM_WORLD);
		remove_files(file, base_dir, myrank);
		//toc_mpi(t_list, 1);
		MPI_Finalize();
	}
	if(bchunk==1){
		//int *ac = (int*)malloc(sizeof(int*));
		//char ***av = (char***)malloc(sizeof(char***));
		//MPI_Init(ac, av);	
		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		//tic_mpi(t_list, 1);
		sint myrank, size;
		MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
		MPI_Comm_size (MPI_COMM_WORLD, &size);
		printf("start removing dup\n");
		remove_duplicate(file, base_dir, myrank);
		MPI_Barrier(MPI_COMM_WORLD);
		if(num_nodes>1){
			//printf("OFFSET: %d\n", OFFSET);
			lint *distrib = (lint*)malloc(sizeof(lint)*num_nodes);
			lint *send_size = (lint*)malloc(sizeof(lint)*num_nodes);
			lint *recv_size = (lint*)malloc(sizeof(lint)*num_nodes);
			comp_distrib_csr_par_bin(file, base_dir, distrib, send_size, recv_size, SEND, myrank);
			MPI_Barrier(MPI_COMM_WORLD);
			int i;
			if(myrank==0) printf("start computing chunk_files_bin_par\n");
			chunk_files_bin_par(file, base_dir, distrib, SEND, myrank);
			chunk_files_bin_par(file, base_dir, distrib, RECV, myrank);
			MPI_Barrier(MPI_COMM_WORLD);
			bin_remove_files(file, base_dir, myrank);
		}
		//toc_mpi(t_list, 1);
		MPI_Finalize();
	}
	if(mpi==1){
		//tic_mpi(t_list, 1);
		run_pagerank_csr_mpi(file, base_dir);
		//toc_mpi(t_list, 1);
	}
	if(delta_stepping == 1){
		//run_sspath_list(file, base_dir, 0, 0.1);
		//run_djkstra(file, base_dir, 0);
		//initialization step
		lint i;
		t_csr *gs = (t_csr*)malloc(sizeof(t_csr));
		scan_csr_idx(gs, file, base_dir, SEND);
		printf("finished scanning graph\n");
		read_graph_csr(gs, file, base_dir, SEND);
		printf("finished scanning graph\n");
		p_list *predecesor = (p_list*)malloc(sizeof(p_list)*gs->v_size);
		for(i=0;i<gs->v_size;i++){
			predecesor[i] = (p_list)malloc(sizeof(t_list));
			creat_list(predecesor[i]);
		}
		printf("finished p_list allocation\n");
		//run_sspath(gs, predecesor, 0, 0.1);
		//run_djkstra(file, base_dir, 0);
		run_sspath_arr(gs, predecesor, 0, 0.4);

	}
	if(adj==1)
	{
		return 0;
	}
	if(bc==1){
		run_bc(file, base_dir);
	}
	if(reduce == 1){
		reduce_graph_for_pr(file, base_dir);
	}
	else if(csr==1)
	{
		run_pagerank_csr_mpi(file, base_dir);
		return 0;
	}
	else if(help==1){
		printf("USAGE: ./PR -<options> [args]\n");
		printf("-a:  using adjacency list\n");
		printf("-b:  setting base directory\n");
		printf("-c:  using csr\n");
		printf("-d:  run delta stepping algorithm\n");
		printf("-e:  run betweeness centrality algorithm\n");
		printf("-f:  setting graph file\n");
		printf("-g:  parallel partition on \n");
		printf("-h:  print help message \n");
		printf("-h:  running with mpi enabled \n");
		printf("-k:  running serial chunking  \n");
		printf("-n:  how many mpi nodes are used  \n");
		printf("-p:  how many omp threads are used  \n");
		printf("-s:  set up the scale \n");
		printf("-y:  parallel chunking on binary files \n");
	}

	//printf("page rank time is %f seconds.\n", t_list->list[0].time);
	return 0;
}

