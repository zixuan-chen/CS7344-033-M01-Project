#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main (int argc, char *argv[]) {
    int count; /* Local prime count */ 
    double elapsed_time; /* Parallel execution time */ 
    int i;
    int j;  
    int id; /* Process ID number */ 
    char *marked; /* Portion of 2,...,'n' */ 
    char *is_prime;
    char *global_marked;
    int n; /* Sieving from 2, ..., 'n' */ 
    int p; /* Number of processes */ 
    int limit; /* Elements in 'marked' */ 

    MPI_Init (&argc, &argv); 
    MPI_Barrier(MPI_COMM_WORLD); 
    elapsed_time = -MPI_Wtime(); 
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    if (argc != 2) 
    { 
        if (!id) 
        printf ("Command line: %s <m>\n", argv[0]);
        MPI_Finalize(); 
        exit (1); 
    }

    n = atoi(argv[1]); 
    limit = sqrt(n);

    marked = (char *) malloc (n+1); 
    global_marked = (char *) malloc (n+1); 
    is_prime = (char *)malloc (limit+1);
    if (marked == NULL || is_prime == NULL) { 
        printf ("Cannot allocate enough memory\n"); 
        MPI_Finalize(); 
        exit (1);
    }
    for (i = 2; i <= limit; i++) is_prime[i] = 1; 
    for (i = 2; i <= n; i++) {marked[i] = 0; global_marked[i] = 0;}

    for (i = 2; i <= (int)(sqrt(limit)); ++i)
    {
        if (is_prime[i])
            for(j = 2*i; j <= limit; j+=i)
                is_prime[j] = 0;
    }

    count = 0;
    for (i = 2; i <= limit; ++i){
        if(is_prime[i]) {
            count++;
            if(count % p == id) {
                printf("process %d catches %d\n", id, i);
                for (j = 2*i; j <= n; j+=i)
                    marked[j] = 1;    
            }
        }
    }
    MPI_Reduce(marked, global_marked, n+1, MPI_CHAR, MPI_LOR, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&count, 1, MPI_INT, id, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime(); 
    if (!id) { 
        count = 0; 
        for (i = 2; i <= n; i++) {
            if (!global_marked[i]) count++;
        }
        printf ("%d primes are less than or equal to %d\n", count, n); 
        printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }
    MPI_Finalize (); 
    return 0;

}