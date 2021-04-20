#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int my_MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    int p, id;
    MPI_Status status;
    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id);
    if (id == root)
    {
        for (int i = 0; i < p; ++i)
            if (i != id)
                MPI_Send(buffer, count, datatype, i, 0, comm);
    }
    else{
        MPI_Recv(buffer, count, datatype, root, 0, comm, &status);
    }
    return 0;
}

int main (int argc, char *argv[]) {
    int *value; /* Local prime count */ 
    double elapsed_time; /* Parallel execution time */ 
    int id; /* Process ID number */ 
    int n; /* Sieving from 2, ..., 'n' */ 
    int p; /* Number of processes */ 
    int i;

    MPI_Init (&argc, &argv); 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    n = atoi(argv[1]); 

    value = (int *) malloc (n*sizeof(int)); 
    if (value == NULL)
    {
        printf ("Cannot allocate enough memory\n"); 
        MPI_Finalize(); 
        exit (1);
    }
    for (i = 0; i <= n; i++) value[i] = id; 
    printf("Using my_MPI_Bcast: process %d is initialized with value = %d\n", id, value[0]);
    elapsed_time = -MPI_Wtime(); 
    my_MPI_Bcast(value, n, MPI_INT, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    printf("Using my_MPI_Bcast: process %d updated with value = %d, time passed: %f\n", id, value[0], elapsed_time);

    
    for (i = 0; i <= n; i++) value[i] = id; 
    printf("Using MPI_Bcast: process %d is re-initialized with value = %d\n", id, value[0]);
    elapsed_time = -MPI_Wtime(); 
    MPI_Bcast(value, n, MPI_INT, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    printf("Using MPI_Bcast: process %d updated with value = %d, time passed: %f\n", id, value[0], elapsed_time);
    MPI_Finalize (); 
    return 0;
}