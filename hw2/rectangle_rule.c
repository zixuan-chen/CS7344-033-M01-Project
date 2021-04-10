#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* This program computes pi using the rectangle rule. */ 
#define INTERVALS 1000000 


#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1) 
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW(id,p,n)) 

int main (int argc, char *argv[])
{ 
    double area; /* Area under curve */ 
    double global_area; /* Global area under curve */
    double ysum; /* Sum of rectangle heights */ 
    double xi; /* Midpoint of interval */ 
    int i; 
    int id; /* Process ID number */ 
    int p; /* Number of processes */ 
    int low_value; /* Lowest value on this proc */ 
    int high_value; /* Highest value on this proc */
    int size; /* Elements in 'marked' */ 
    double elapsed_time; /* Parallel execution time */ 

    MPI_Init (&argc, &argv); 
    MPI_Barrier(MPI_COMM_WORLD); 
    elapsed_time = -MPI_Wtime(); 
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 

    low_value = BLOCK_LOW(id,p,INTERVALS); 
    high_value = BLOCK_HIGH(id,p,INTERVALS); 
    size = BLOCK_SIZE(id,p,INTERVALS); 
    ysum = 0.0; 
    for (i = low_value; i <= high_value; i++) { 
        xi = (1.0/INTERVALS)*(i+0.5); 
        ysum += 4.0/(1.0+xi*xi); 
    } 
    area = ysum * (1.0 / INTERVALS); 

    // printf("process No.%d, low_value=%d, high_value=%d, size=%d, area = %10.6f\n", 
    //         id, low_value, high_value, size, area);
    MPI_Reduce (&area, &global_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime(); 
    if (!id)
    {
        printf ("Area is %13.11f\n", global_area); 
        printf ("Total elapsed time: %10.6f\n", elapsed_time);
    }
    MPI_Finalize();
    return 0; 
} 


