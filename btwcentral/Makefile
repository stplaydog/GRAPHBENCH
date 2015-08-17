include ../Makefile.common

csr_bc.a: csr_bc.o
	$(AR) rcs ../lib/csr_bc.a csr_bc.o

csr_bc.o: csr_bc.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr_bc.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/csr_bc.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/csr_bc.a
