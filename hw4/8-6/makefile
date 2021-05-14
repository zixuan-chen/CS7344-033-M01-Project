
CC = mpicc
CFLAGS = -Wall -g3

executables = exer08_06


all : $(executables)


# executable construction ----------------------------------

$(executables) : $(executables).o MyMPI.o
	$(CC) $(CFLAGS) $(executables).o MyMPI.o -lm -o $(executables)


# object file construction ---------------------------------


$(executables).o : $(executables).c MyMPI.h
	$(CC) $(CFLAGS) -c $(executables).c


MyMPI.o : MyMPI.c MyMPI.h
	$(CC) $(CFLAGS) -c MyMPI.c



# file cleanup ---------------------------------------------

.PHONY : clean clean_o clean_execut

clean :
	rm *.o $(executables) -f

clean_o :
	rm *.o -f

clean_execut :
	rm $(executables) -f

