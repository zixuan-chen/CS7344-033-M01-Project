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
    int index; /* Index of current prime */ 
    int low_value; /* Lowest value on this proc */ 
    char *marked; /* Portion of 2,...,'n' */ 
    int n; /* Sieving from 2, ..., 'n' */ 
    int p; /* Number of processes */ 
    int proc0_size; /* Size of proc 0's subarray */ 
    int prime; /* Current prime */ 
    int size; /* Elements in 'marked' */ 

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
    low_value = 2 + BLOCK_LOW(id,p,n-1); 
    high_value = 2 + BLOCK_HIGH(id,p,n-1); 
    size = BLOCK_SIZE(id, p, n-1); 

    /* Bail out if all the primes used for sieving are not all held by process 0 */ 
    proc0_size = (n-1)/p; 
    if ((2 + proc0_size) < (int) sqrt((double) n)) { 
        if (!id) 
            printf ("Too many processes\n"); 
        MPI_Finalize(); 
        exit (1);
    }

    marked = (char *) malloc (size + 2); 
    if (marked == NULL) { 
        printf ("Cannot allocate enough memory\n"); 
        MPI_Finalize(); 
        exit (1);
    }
    for (i = 0; i < size + 2; i++) marked[i] = 0; 
    if (!id) index = 0; 
    prime = 2; 
    do 
    { 
        if (prime * prime > low_value) 
            first = prime * prime - low_value; 
        else { 
            if (!(low_value % prime)) first = 0; 
            else first = prime - (low_value % prime);
        }
        /* check primes */
        for (i = first; i < size + 2; i += prime) marked[i] = 1; 
        if (!id) { 
            while (marked[++index]); 
            prime = index + 2;
        }
        MPI_Bcast (&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);

    } while (prime * prime <= n);
    count = 0; 
    consecutive_count = 0;
    for (i = 0; i < size; i++) 
    {
        if (!marked[i]) count++;
        if (!marked[i] && !marked[i+2]) consecutive_count++;
    }
        

    MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce (&consecutive_count, &global_consecutive_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime(); 
    if (!id) { 
        printf ("%d primes are less than or equal to %d\n", global_count, n); 
        printf ("among which, %d primes are consecutive\n", global_consecutive_count); 
        printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }
    MPI_Finalize (); 
    return 0;

}