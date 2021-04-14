#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1) 
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW(id,p,n)) 
#define BLOCK_OwNER(index,p,n) (((p)*((index)+1)-1)/(n)) 

int main (int argc, char *argv[]) {
    int count; /* Local prime count */ 
    int consecutive_count; /* Consecutive local prime count */
    double elapsed_time; /* Parallel execution time */ 
    int first; /* Index of first multiple */ 
    int global_count; /* Global prime count */ 
    int global_consecutive_count; /* Global consecutive prime count */ 
    int high_value; /* Highest value on this proc */
    int i;  
    int id; /* Process ID number */ 
    int low_value; /* Lowest value on this proc */ 
    char *marked; /* Portion of 2,...,'n' */ 
    int n; /* Sieving from 2, ..., 'n' */ 
    int p; /* Number of processes */ 
    int proc0_size; /* limit of proc 0's subarray */ 
    int prime; /* Current prime */ 
    int limit; /* Elements in 'marked' */ 
    int fire;

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
    if (marked == NULL) { 
        printf ("Cannot allocate enough memory\n"); 
        MPI_Finalize(); 
        exit (1);
    }
    for (i = 2; i <= n; i++) marked[i] = 0; 
    prime = 2;
    fire = 0;
    while(prime <= limit)
    {
        if(!id)
        {
            while(prime <= limit && marked[prime]) prime++;
            fire = (fire+1) % p; // decide which process to fire
            MPI_Bcast (&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast (&fire, 1, MPI_INT, 0, MPI_COMM_WORLD);    
            printf("process %d: fire=%d, prime=%d\n", id, fire, prime);
        }
        
        if (id == fire)
        {
            printf("process %d firing: sieving on prime %d\n", id, prime);
            for(i = prime; i <= n; i+=prime)
                marked[i] = 1;
            // prime++;
            // MPI_Bcast (&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
            /*third step: OR_reduce*/
            if (id)
                MPI_Reduce(&marked, &marked, n+1, MPI_CHAR, MPI_LOR, 0, MPI_COMM_WORLD);
        }
    }
    
    elapsed_time += MPI_Wtime(); 
    if (!id) { 
        count = 0; 
        for (i = 2; i <= n; i++) 
        {
            if (!marked[i]) count++;
        }
        printf ("%d primes are less than or equal to %d\n", count, n); 
        printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }
    MPI_Finalize (); 
    return 0;

}