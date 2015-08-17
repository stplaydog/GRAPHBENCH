include ../Makefile.common

sssp.a: csr_sspath_array.o csr_sspath.o
	$(AR) rcs ../lib/sssp.a csr_sspath_array.o csr_sspath.o

csr_sspath_array.o: csr_sspath_array.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_sspath_array.c 
csr_sspath.o: csr_sspath.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_sspath.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/sssp.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/sssp.a
