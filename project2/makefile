OBJS = monte_carlo pagerank quicksort
cc = g++
CC = mpiCC
CFLAGS = -Wall

.PHONY: all, clean


all: $(OBJS)
	@echo Done!


clean:
	rm -f $(OBJS)

clear:
	rm -f $(OBJS)


monte_carlo: monte_carlo.cpp
	$(cc) $< -o $@ -fopenmp

pagerank: pagerank.cpp
	$(CC) $< -o $@ -fopenmp

quicksort: quicksort.cpp
	$(CC) $< -o $@ -fopenmp