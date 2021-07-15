OBJS = my_MPI_Allgather matmul conv wordcount_small wordcount_big
cc = mpicc
CC = mpiCC
CFLAGS = -Wall

.PHONY: all, clean


all: $(OBJS)
	@echo Done!


clean:
	rm -f $(OBJS)

clear:
	rm -f $(OBJS)

# matrix: gen_matrix.py 
# 	python3 gen_matrix.py

my_MPI_Allgather: my_MPI_Allgather.c
	$(cc) $< -o $@

matmul: matmul.c MyMPI.c MyMPI.h 
	$(cc) -D PRINT_FLAG=0 MyMPI.c $< -o $@

conv: conv.c  MyMPI.c MyMPI.h
	$(cc) -D PRINT_FLAG=0 MyMPI.c  $< -o $@

wordcount_small: wordcount_small.cpp wordcount.h MyMPI.c MyMPI.h
	$(CC) $< -o $@

wordcount_big: wordcount_big.cpp wordcount.h MyMPI.c MyMPI.h
	$(CC) $< -o $@