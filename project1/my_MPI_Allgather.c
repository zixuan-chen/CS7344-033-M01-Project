#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int my_MPI_Allgather(const void* sendbuf, int sendcount, 
                    MPI_Datatype sendtype, void* recvbuf, int recvcount,
                    MPI_Datatype recvtype, MPI_Comm comm){
    int p, id;
    MPI_Status status;
    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id);
    int i, j;
    int sendsize, recvsize;
    void *recvptr;
    MPI_Type_size(recvtype, &recvsize);
    MPI_Type_size(sendtype, &sendsize);
    
    for(i = 0; i < p; ++i){
        if (i  != id)
            MPI_Send(sendbuf, sendcount, sendtype, i, 0, comm);
        else{
            memcpy(recvbuf+i*recvcount*recvsize, sendbuf, sendcount*sendsize);
        }
    }
    
    recvptr = recvbuf;
    for(i = 0; i < p; ++i){
        if (i != id)
            MPI_Recv(recvptr, recvcount, recvtype, i, 0, comm, &status);
        recvptr = recvptr + recvsize*recvcount;
    }
        
    return 0;
}

int main (int argc, char *argv[]) {
    int *send_buff; /* Local prime count */ 
    int *recv_buff; 
    double elapsed_time; /* Parallel execution time */ 
    int id; /* Process ID number */ 
    int n; 
    int p; /* Number of processes */ 
    int i;

    MPI_Init (&argc, &argv); 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    n = atoi(argv[1]); 
    int size;
    
    send_buff = (int *) malloc (n*sizeof(int)); 
    recv_buff = (int *) malloc (n*sizeof(int));
    if (send_buff == NULL || recv_buff == NULL)
    {
        printf ("Cannot allocate enough memory\n"); 
        MPI_Finalize(); 
        exit (1);
    }
    // printf("Using my_MPI_Bcast: process %d is initialized with value = %d\n", id, send_buff[0]);
    // elapsed_time = -MPI_Wtime(); 
    // my_MPI_Bcast(send_buff, n, MPI_INT, 0, MPI_COMM_WORLD);
    // elapsed_time += MPI_Wtime();
    // printf("Using my_MPI_Bcast: process %d updated with value = %d, time passed: %f\n", id, send_buff[0], elapsed_time);
    elapsed_time = 0;
    if(!id)
    printf("Using my_MPI_Allgather: \n");
    elapsed_time = -MPI_Wtime(); 
    my_MPI_Allgather(send_buff, n/p, MPI_INT, recv_buff, n/p, MPI_INT, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(!id)printf("Total elapsed time: %10.6f\n", elapsed_time);

    elapsed_time = 0;
    if(!id)printf("Using MPI_Allgather: \n");
    elapsed_time = -MPI_Wtime(); 
    MPI_Allgather(send_buff, n/p, MPI_INT, recv_buff, n/p, MPI_INT, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if(!id)printf("Total elapsed time: %10.6f\n", elapsed_time);

    
    MPI_Finalize (); 
    return 0;
}