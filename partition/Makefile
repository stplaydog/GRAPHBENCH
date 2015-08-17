include ../Makefile.common

partition.a: csr_partition_mpi.o csr_partition.o csr_partition_bin_mpi.o
	$(AR) rcs ../lib/partition.a csr_partition_mpi.o csr_partition.o csr_partition_bin_mpi.o

csr_partition_bin_mpi.o: csr_partition_bin_mpi.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_partition_bin_mpi.c 
csr_partition.o: csr_partition.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_partition.c 
csr_partition_mpi.o: csr_partition_mpi.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_partition_mpi.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/partition.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/partition.a
