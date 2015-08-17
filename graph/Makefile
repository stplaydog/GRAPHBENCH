include ../Makefile.common

csr.a: csr.o
	$(AR) rcs ../lib/csr.a csr.o

csr.o: csr.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c csr.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/csr.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/csr.a
