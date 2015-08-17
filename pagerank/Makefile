include ../Makefile.common

pagerank.a: pagerank_mpi.o preprocess.o
	$(AR) rcs ../lib/pagerank.a pagerank_mpi.o preprocess.o

pagerank_mpi.o: pagerank_mpi.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c pagerank_mpi.c 
preprocess.o: preprocess.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c preprocess.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/pagerank.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/pagerank.a
