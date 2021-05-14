#include <stdio.h>
#include <mpi.h>
#include "MyMPI.h"

typedef double dtype;
#define mpitype MPI_DOUBLE


int main(int argc, char *argv[]) 
{
    dtype **a;              /* The first factor, a matrix */ 
    dtype  *b;              /* The second factor, a vector */ 
    dtype  *c;              /* The product, a vector */ 
    dtype  *c_part_out;     /* Partial sums, sent */ 
    dtype  *c_part_in;      /* Partial sums, received */ 
    int    *cnt_out;        /* Elements sent to each proc */ 
    int    *cnt_in;         /* Elements received per proc */ 
    int    *disp_out;       /* Indices of sent elements */ 
    int    *disp_in;        /* Indices of received elements */ 
    int     i, j;           /* Loop indices */ 
    int     id;             /* Process ID number */ 
    int     local_low;      
    int     local_els;      /* Cols of 'a' and elements of 'b' held by this process */ 
    int     m;              /* Rows in the matrix */ 
    int     n;              /* Columns in the matrix */ 
    int     nprime;         /* Size of the vector */ 
    int     p;              /* Number of processes */ 
    dtype  *storage;        /* This process's portion of 'a' */

    MPI_Init (&argc, &argv); 
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    m = 0; n = 0;
    read_col_striped_matrix (argv[1], (void ***) &a, (void **) &storage, 
                                mpitype, &m, &n, MPI_COMM_WORLD); 
    print_col_striped_matrix ((void **) a, mpitype, m, n, MPI_COMM_WORLD); 
    read_replicated_vector (argv[2], (void **) &b, mpitype, &nprime, MPI_COMM_WORLD); 
    print_replicated_vector ((void *) b, mpitype, nprime, MPI_COMM_WORLD); 

    /*  
        Each process multiplies its columns of 'a' and vector 'b', 
        resulting in a 
        partial sum of product 'c'. 
    */ 

    c_part_out = (dtype *) my_malloc (id, m * sizeof(dtype)); 
    
    for (i = 0; i < m; i++) { 
        c_part_out[i] = 0.0; 
        for (j = 0; j < BLOCK_SIZE(id,p,n); j++) 
            c_part_out[i] += a[i][j] * b[BLOCK_LOW(id,p,n) + j]; 
    } 

    local_els = BLOCK_SIZE(id,p,m);
    create_mixed_xfer_arrays (id, p, m, &cnt_out, &disp_out); 
    create_uniform_xfer_arrays (id, p, m, &cnt_in, &disp_in); 
    c_part_in = (dtype*) my_malloc (id, p*local_els*sizeof(dtype)); 
    MPI_Alltoallv (c_part_out, cnt_out, disp_out, mpitype, c_part_in, 
                    cnt_in, disp_in, mpitype, MPI_COMM_WORLD); 
    c = (dtype*) my_malloc (id, local_els * sizeof(dtype)); 
    
    for (i = 0; i < local_els; i++) { 
        c[i] = 0.0; 
        for (j = 0; j < p; j++) 
            c[i] += c_part_in[i + j*local_els]; 
    } 
    print_block_vector ((void *) c, mpitype, m, MPI_COMM_WORLD); 
    MPI_Finalize(); 
    return 0; 
}