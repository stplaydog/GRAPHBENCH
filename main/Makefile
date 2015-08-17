include ../Makefile.common

SNY: main.o ../lib/csr_bc.a ../lib/csr_pagerank.a ../lib/csr_sssp.a ../lib/partition.a ../lib/utils.a ../lib/structs.a
	$(CC) $(FLAGS) $(LIB) $(INC) -o SNY main.o ../lib/csr_bc.a ../lib/csr_pagerank.a ../lib/csr_sssp.a ../lib/partition.a ../lib/utils.a ../lib/structs.a$(LINK)

main.o: main.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c main.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/main.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/main.a
