#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define FILE_NAME "input.txt"
#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)) 
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW(id,p,n)) 
#define BLOCK_OwNER(index,p,n) (((p)*((index)+1)-1)/(n)) 


int count_grid(const int * grid, int i, int m, int n)
{
    int count;
    int x = i/m, y = i%m;
    int xx, yy;
    int ii, jj;
    
    count = 0;
    for(ii = -1; ii <= 1; ++ii)
        for(jj = -1; jj <=1; ++jj){
            xx = x + ii;
            yy = y + jj;
            if(xx>=0 && xx<m && yy>=0 && yy<n)
                count += grid[xx*n+yy];
        }
    count -= grid[i];
    //printf("pos = %d, x = %d, y = %d, count = %d\n", i, x, y, count);
    return count;
}

int main (int argc, char *argv[]) {
    double elapsed_time; /* Parallel execution time */ 
    int i, j, m, n, p, k;  
    int id; /* Process ID number */ 
    int* cell_state;
    int* count;
    MPI_Status status;
    int x, y;
    int low, high, size;
    int c;
    

    FILE* fp = fopen(FILE_NAME, "r");
    if (fp == NULL || fscanf(fp, "%d %d", &m, &n) == EOF)
        printf("ERROR: read from file %s.\n", FILE_NAME);
    if (m > 0 && n >0){
        cell_state = (int*)malloc(m*n*sizeof(int));
        count = (int*)malloc(m*n*sizeof(int));
    }
    else
        printf("ERROR: m = %d, n = %d\n", m, n);
    for(x = 0; x < m; ++x)
        for(y = 0; y < n; ++y)
            fscanf(fp, "%d", &cell_state[x*n+y]);
    fclose(fp);

    MPI_Init (&argc, &argv); 
    MPI_Barrier(MPI_COMM_WORLD); 
    elapsed_time = -MPI_Wtime(); 
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    if (p > m*n/2)
    {
        if(!id)
            printf("too many processes\n");
        MPI_Finalize(); 
        exit (1); 
    }
    if (argc != 3) 
    { 
        if (!id) 
        printf ("Command line: %s <m>\n", argv[0]);
        MPI_Finalize(); 
        exit (1); 
    }
    j = atoi(argv[1]);
    k = atoi(argv[2]);

    low = BLOCK_LOW(id, p, m*n);
    high = BLOCK_HIGH(id, p, m*n);
    size = BLOCK_SIZE(id, p, m*n);
    // printf("low = %d, high = %d, size = %d\n", low, high, size);

    for(i = 0; i < j; ++i)
    {   
        MPI_Bcast(cell_state, m*n, MPI_INT, 0, MPI_COMM_WORLD);
        if (!id && (i+1) % k == 0){
            printf("cell state at iteration %d:\n", i+1);
            for(x = 0; x < m; ++x){
                for(y = 0; y < n; ++y)
                    printf("%d ", cell_state[x*n+y]);
                printf("\n");
            }   
        }
        for (c = low; c < high; ++c)
            count[c] = count_grid(cell_state, c, m, n);
        
        for (c = low; c < high; ++c){
            if(!cell_state[c] && count[c] == 3)
                cell_state[c] = 1;
            else if (cell_state[c] && (count[c] < 2 || count[c] > 3))
                cell_state[c] = 0;
        }
        if(id)
            MPI_Send(cell_state+low, size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        else
            {
            for (int pid = 1; pid < p; ++pid)
                MPI_Recv(cell_state+BLOCK_LOW(pid, p, m*n), 
                        BLOCK_SIZE(pid, p, m*n), 
                        MPI_INT, 
                        pid, 
                        0, 
                        MPI_COMM_WORLD, 
                        &status);
            }
    }
    
    MPI_Finalize (); 
    return 0;

}