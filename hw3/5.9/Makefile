CC = mpicc
CFLAGS = -Wall -g3

executables = eratosthenes

all : $(executables)


# executable construction ----------------------------------


$(executables) : $(executables).o
	$(CC) $(CFLAGS) $(executables).o -o $(executables) -lm



# object file construction ---------------------------------

$(executables).o : $(executables).c
	$(CC) $(CFLAGS) -c $(executables).c


# file cleanup ---------------------------------------------

.PHONY : clean clean_o clean_execut

clean :
	rm *.o $(executables) -f

clean_o :
	rm *.o -f

clean_execut :
	rm $(executables) -f