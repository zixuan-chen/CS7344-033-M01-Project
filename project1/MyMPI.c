#include <stdio.h> 
#include <stdlib.h> 
#include <mpi.h> 
#include "MyMPI.h" 
/***************** MISCELLANEOUS FUNCTIONS ********** ***** **/ 

/* 
• Given MPI_Datatype 't', function 'get_size' returns the 
• size of a single datum of that data type. 
*/ 

int get_size (MPI_Datatype t) 
{
    if (t == MPI_BYTE) return sizeof(char); 
    if (t == MPI_DOUBLE) return sizeof(double); 
    if (t == MPI_FLOAT) return sizeof(float); 
    if (t == MPI_INT) return sizeof(int); 
    else{
        printf ("Error: Unrecognized argument to 'get_size'\n"); 
        fflush (stdout); 
        MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR); 
    }
}
/* 
• Function 'my_malloc' is called when a process wants 
• to allocate some space from the heap. If the memory 
* allocation fails, the process prints an error message 
* and then aborts execution of the program. 
*/ 

void *my_malloc (
    int id,         /* IN - Process rank */ 
    int bytes)      /* IN - Bytes to allocate */ 
{   
    void *buffer; 
    if ((buffer = malloc (bytes)) == NULL) 
    { 
        printf ("Error: Malloc failed for process %d\n", id); 
        fflush (stdout); 
        MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR); 
    }
    return buffer; 
} 


/* 
• Function 'terminate' is called when the program should 
• not continue execution, due to an error condition that all of the processes are aware of. 
* Process 0 prints the error message passed as an argument to the function. 
* All processes must invoke this function together! */ 

void terminate (
    int id, /* IN - Process rank */ 
    char *error_message) /* IN - Message to print */ 
{
    if (!id) printf ("Error: %s\n", error_message); 
    fflush (stdout); 
    MPI_Finalize(); 
    exit (-1); 
}

/***************DATA DISTRIBUTION FUNCTIONS******************************/

/* 
* This function creates the count and displacement arrays 
* needed by scatter and gather functions, when the number 
* of elements send/received to/from other processes 
* varies. 
*/

void create_mixed_xfer_arrays (
    int     id,     /* IN - Process rank */ 
    int     p,      /* IN - Number of processes */ 
    int     n,      /* IN - Total number of elements */ 
    int   **count,  /* OUT - Array of counts */ 
    int   **disp)   /* OUT - Array of displacements */ 
{
    int i; 
    *count = my_malloc (id, p * sizeof(int)); 
    *disp = my_malloc (id, p * sizeof(int)); 
    (*count)[0] = BLOCK_SIZE(0,p,n); 
    (*disp)[0] = 0; 
    for (i = 1; i < p; i++) {
        (*disp)[i]=(*disp)[i-1] + (*count)[i-1]; 
        (*count)[i] = BLOCK_SIZE(i,p,n); 
    }
}

/* 
* This function creates the count and displacement arrays 
• needed in an all-to-all exchange, when a process gets 
* the same number of elements from every other process. 
*/
void create_uniform_xfer_arrays (
    int     id,     /* IN - Process rank */ 
    int     p,      /* IN - Number of processes */ 
    int     n,      /* IN - Number of elements */ 
    int   **count,  /* OUT - Array of counts */ 
    int   **disp)   /* OUT - Array of displacements */ 
{
    int i; 
    *count = my_malloc (id, p * sizeof(int)); 
    *disp = my_malloc (id, p * sizeof(int)); 
    (*count)[0] = BLOCK_SIZE(id,p,n); 
    (*disp)[0] = 0; 
    for (i = 1; i < p; i++) { 
        (*disp)[i]=(*disp)[i-1] + (*count)[i-1]; 
        (*count)[i] = BLOCK_SIZE(id,p,n); 
    }
}
/* 
* This function is used to transform a vector from a 
* block distribution to a replicated distribution within a 
* * communicator. 
*/

void replicate_block_vector (
    void *ablock, /* IN - Block-distributed vector */ 
    int n, /* IN - Elements in vector */ 
    void *arep, /* OUT - Replicated vector */ 
    MPI_Datatype dtype, /* IN - Element type */ 
    MPI_Comm comm) /* IN - Communicator */ 
{
    int *cnt; /* Elements contributed by each process */ 
    int *disp; /* Displacement in concatenated array */ 
    int id; /* Process id */ 
    int p; /* Processes in communicator */ 
    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id); 
    create_mixed_xfer_arrays (id, p, n, &cnt, &disp); 
    MPI_Allgatherv (ablock, cnt[id], dtype, arep, cnt, disp, dtype, comm);
    free (cnt); 
    free (disp);  
}   
/* 
• Function 'read_checkerboard_matrix' reads a matrix from 
• a file. The first two elements of the file are integers 
• whose values are the dimensions of the matrix ('m' rows • and 'n' columns). 
What follows are 'm'*'n' values • representing the matrix elements stored in row-major 
• order. This function allocates blocks of the matrix to the MPI processes. 
* The number of processes must be a square number. 
*/ 
void read_checkerboard_matrix(
    char *s,                /* IN - File name */
    void ***subs,           /* OUT - 2D array */
    void **storage,         /* OUT - Array elements */
    MPI_Datatype dtype,     /* IN - Element type */
    int *m,                 /* OUT - Array rows */
    int *n,                 /* OUT - Array cols */
    MPI_Comm grid_comm)     /* IN - Communicator */ 
{
    void *buffer;           /* File buffer */ 
    int coords[2];          /* Coords of proc receiving next row of matrix */ 
    int datum_size;         /* Bytes per element */ 
    int dest_id;            /* Rank of receiving proc */ 
    int grid_coord[2];      /* Process coords */ 
    int grid_id;            /* Process rank */ 
    int grid_period[2];     /* Wraparound */ 
    int grid_size[2];       /* Dimensions of grid */ 
    int i, j, k; 
    FILE *infileptr;        /* Input file pointer */ 
    void *laddr;            /* Used when proc 0 gets row */ 
    int local_cols;         /* Matrix cols on this proc */ 
    int local_rows;         /* Matrix rows on this proc */ 
    void **lptr;            /* Pointer into 'subs' */ 
    int p;                  /* Number of processes */ 
    void *raddr;            /* Address of first element to send */ 
    void *rptr;             /* Pointer into 'storage' */ 
    MPI_Status status;      /* Results of read */ 
    MPI_Comm_rank (grid_comm, &grid_id); 
    MPI_Comm_size (grid_comm, &p); 
    datum_size = get_size (dtype); 
    /* Process 0 opens file, gets number of rows and number of cols, 
    and broadcasts this information to the other processes. */ 
    if (grid_id == 0) {
        infileptr = fopen (s, "r"); 
        if (infileptr == NULL) *m = 0; 
        else {
            fscanf(infileptr, "%d", m);
            fscanf(infileptr, "%d", n); 
        }
    }
    MPI_Bcast (m, 1, MPI_INT, 0, grid_comm); 
    if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR); 
    MPI_Bcast (n, 1, MPI_INT, 0, grid_comm); 
    /* Each process determines the size of 
    the submatrix it is responsible for. */ 
    MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coord); 
    local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m); 
    local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n); 
    /* Dynamically allocate two-dimensional matrix 'subs' */ 
    *storage = my_malloc (grid_id, local_rows * local_cols * datum_size); 
    *subs = (void **) my_malloc (grid_id,local_rows*PTR_SIZE); 
    lptr = (void *) *subs; 
    rptr = (void *) *storage; 
    for (i = 0; i < local_rows; i++) {
        *(lptr++) = (void *) rptr; 
        rptr += local_cols * datum_size; 
    }
    /* Grid process 0 reads in the matrix one row at a time 
    and distributes each row among the MPI processes. */ 
    if (grid_id == 0) 
        buffer = my_malloc (grid_id, *n * datum_size); 
    /* For each row of processes in the process grid... */ 
    for (i = 0; i < grid_size[0]; i++) { 
        coords[0] = i; 
        /* For each matrix row controlled by this proc row...*/ 
        for (j = 0; j < BLOCK_SIZE(i,grid_size[0],*m); j++) { 
            /* Read in a row of the matrix */ 
            if (grid_id == 0) { 
                for (k = 0; k < *n; ++ k) {
                    fscanf(infileptr, SCANF_FORMAT, (D_TYPE *)(buffer + k*datum_size));
                }
            }

            /* Distribute it among processes in the grid row */ 
            for (k = 0; k < grid_size[1]; k++) { 
                coords[1] = k; 
                /* Find address of first element to send */ 
                raddr = buffer + BLOCK_LOW(k,grid_size[1],*n) * datum_size; 
                /* Determine the grid ID of the process getting the subrow */ 
                MPI_Cart_rank (grid_comm, coords, &dest_id); 

                /* Process 0 is responsible for sending...*/ 
                if (grid_id == 0) { 
                    /* It is sending (copying) to itself */ 
                    if (dest_id == 0) { 
                        laddr = (*subs)[j]; 
                        memcpy (laddr, raddr, local_cols * datum_size);
                    /* It is sending to another process */ 
                    } else { 
                        MPI_Send (raddr, 
                                BLOCK_SIZE(k,grid_size[1],*n), 
                                dtype, dest_id, 0, grid_comm); 
                    }
                /* Process 'dest_id' is responsible for receiving... */ 
                } else if (grid_id == dest_id) {
                    MPI_Recv ((*subs)[j], local_cols, dtype, 0, 0, grid_comm, &status); 
                }
            }
        }
    }
    if(grid_id== 0) free(buffer);
}     
    
void read_simple_matrix(
    char * s,               /* IN - Input file name */
    void ***subs,           /* OUT - 2d matrix pointer */
    void **storage,         /* OUT - 1d matrix storage */
    MPI_Datatype dtype,     /* IN - Element data type */
    int *m,                 /* OUT - Number of rows */
    int *n,                 /* OUT - Number of columns */
    MPI_Comm comm           /* IN - Communicator */
) {
    int id;
    int p;
    int i, j;
    int datum_size;
    FILE *infile_ptr;
    
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);
    MPI_Type_size(dtype, &datum_size);

    if (id == 0) {
        infile_ptr = fopen(s, "r");
        if (infile_ptr == NULL) {
            printf("[error] Open file error\n");
            MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
        }

        fscanf(infile_ptr, "%d", m);
        fscanf(infile_ptr, "%d", n);
    }

    MPI_Bcast(m, 1, MPI_INT, 0, comm);

    if (!(*m)) MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);

    MPI_Bcast(n, 1, MPI_INT, 0, comm);

    *storage = my_malloc(id, (*m) * (*n) * datum_size);
    *subs = (void **)my_malloc(id, (*m) * PTR_SIZE);
    for (i = 0; i < (*m); ++ i) {
        (*subs)[i] = *storage + i * (*n) * datum_size;
    }

    if (id == 0) {
        for (i = 0; i < (*m) * (*n); ++ i) {
            fscanf(infile_ptr, SCANF_FORMAT, (D_TYPE *)((*storage)+i*datum_size));
        }
    }

    MPI_Bcast(*storage, (*m)*(*n), dtype, 0, comm);

    if (id == 0)
        fclose (infile_ptr);
} 
/* Function 'read_col_striped_matrix' reads a matrix from a file.
 The first two elements of the file are integers whose values
  are the dimensions of the matrix ('m' and 'n' columns). 
  What follows are 'm'*'n' values representing the matrix elements 
  stored in row-major order. This function allocates blocks 
  of columns of the matrix to the MPI processes. 
*/
void read_col_striped_matrix(
    char *s,            /* IN -File name */
    void ***subs,       /* OUT - 2-D array */
    void **storage,     /*OUT - Array elements */ 
    MPI_Datatype dtype, /* IN -Element type */
    int *m,             /* 'I OUT - Rows */
    int *n,             /*  OUT - Cols */
    MPI_Comm comm)      /*  IN - Communicator */ 
{
    void *buffer;       /* File buffer */ 
    int datum_size;     /* Size of matrix element */
    int i, j; 
    int id;             /* Process rank */ 
    FILE *infileptr;    /* Input file ptr */ 
    int local_cols;     /* Cols on this process */ 
    void **lptr;        /* Pointer into 'subs' */ 
    void *rptr;         /* Pointer into 'storage' */ 
    int p;              /* Number of processes */ 
    int *send_count;    /* Each proc's count */ 
    int *send_disp;     /* Each proc's displacement */ 

    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id); 
    datum_size = get_size (dtype); 
    /* Process p-1 opens file, gets number of rows and cols, 
    and broadcasts this info to other procs. */ 
    // printf("input matrix file name = %s\n", s);
    if (id == (p-1)) {
        // printf("get into process: %d\n", p-1);
        infileptr = fopen (s, "r"); 
        if (infileptr == NULL) 
        {   
            *m = 0; 
            printf("failed to open file: %s\n", s);
        }
        else 
        { 
            fscanf(infileptr, "%d %d", m, n);
            printf("m = %d, n = %d \n", *m, *n);
        }
    }
    
    MPI_Bcast (m, 1, MPI_INT, p-1, comm); 
    if (!(*m)) 
        MPI_Abort (comm, OPEN_FILE_ERROR); 
    MPI_Bcast (n, 1, MPI_INT, p-1, comm); 
    local_cols = BLOCK_SIZE(id,p,*n); 
    
    /* Dynamically allocate two-dimensional matrix 'subs' */ 
    
    *storage = my_malloc (id, *m * local_cols * datum_size); 
    *subs = (void **) my_malloc (id, *m * PTR_SIZE); 
    lptr = (void *) *subs; 

    rptr = (void *) *storage; 
    for (i = 0; i < *m; i++) { 
        *(lptr++) = (void *) rptr; 
        rptr += local_cols * datum_size; 
    }
    /* 
    Process p-1 reads blocks of rows from file and sends
    each block to the correct destination process. The 
    last block it keeps. 
    */ 
    if (id == (p-1)) 
        buffer = my_malloc (id, *n * datum_size); 
    create_mixed_xfer_arrays (id, p, *n, &send_count, &send_disp); 
    for (i = 0; i < *m; i++) { 
        if (id == (p-1)) {
            for (j = 0; j < *n; j++)
                fscanf(infileptr, "%lf", (double*)buffer+j);
        }
        MPI_Scatterv (buffer, send_count, send_disp, dtype, 
                    (*storage)+i*local_cols*datum_size, local_cols, 
                    dtype, p-1, comm); 
    } 
    free (send_count); 
    free (send_disp); 
    if (id == (p-1)) free (buffer); 
}

/* 
* Open a file containing a vector, read its contents, 
* and distribute the elements by block among the 
* processes in a communicator. 
*/


void read_row_striped_matrix (
    char  *s,           /* IN -  File name */
    void ***subs,       /* OUT - 2D submatrix indices */
    void **storage,     /* OUT - Submatrix stored here */  
    MPI_Datatype dtype, /* OUT - Matrix element type */
    int *m,             /* OUT - Matrix rows */
    int *n,             /* OUT - - Matrix cols */ 
    MPI_Comm comm)      /* IN - Communicator */ 
{ 
    int         datum_size;  /* Size of matrix element */ 
    int         i;
    int         id;         /* Process rank */ 
    FILE        *infileptr; /* Input file pointer */ 
    int         local_rows; /* Rows on this proc */ 
    void        **lptr;     /* Pointer into 'subs' */ 
    int         p;          /* Number of processes */ 
    void        *rptr;      /* Pointer into 'storage' */ 
    MPI_Status  status;     /* Result of receive */ 
    int         x;          /* Result of read */ 

    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id); 
    datum_size = get_size (dtype); 

    /* Process p-1 opens file, reads size of matrix, 
    and broadcasts matrix dimensions to other procs */ 
    if (id == (p-1)) { 
        infileptr = fopen (s, "r"); 
        if (infileptr == NULL) *m = 0; 
        else { 
            fread (m, sizeof(int), 1, infileptr); 
            fread (n, sizeof(int), 1, infileptr); 
        } 
    }
    MPI_Bcast (m, 1, MPI_INT, p-1, comm); 
    if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR); 
    MPI_Bcast (n, 1, MPI_INT, p-1, comm); 
    local_rows = BLOCK_SIZE(id,p,*m); 
    /* Dynamically allocate matrix. Allow double subscripting 
        through 'a'. */ 
    *storage = (void *) my_malloc (id, local_rows * *n * datum_size); 
    *subs = (void **) my_malloc (id, local_rows * PTR_SIZE); 

    lptr = (void *) &(*subs[0]); 
    rptr = (void *) *storage; 
    for (i = 0; i < local_rows; i++) {
        *(lptr++)= (void *) rptr; 
        rptr += *n * datum_size; 
    }

    /* Process p-1 reads blocks of rows from file 
        and sends each block to the correct destination process. 
        The last block it keeps. 
    */ 

    if (id == (p-1)) {
        for (i = 0; i < p-1; i++) { 
            x = fread (*storage, datum_size, BLOCK_SIZE(i,p,*m) * *n, infileptr); 
            MPI_Send (*storage, BLOCK_SIZE(i,p,*m) * *n, dtype, i,  DATA_MSG, comm); 

        } 
        x = fread (*storage, datum_size, local_rows * *n, infileptr); 
        fclose (infileptr); 
    } else 
        MPI_Recv (*storage, local_rows * *n, dtype, p-1, DATA_MSG, comm, &status); 
}
/* 
Open a file containing a vector, read its contents, 
and distribute the elements by block among the processes
in a communicator. 
*/
void read_block_vector(
    char *s,             /* IN - File name */
    void **v,            /* OUT - Subvector */
    MPI_Datatype dtype,  /* IN - Element type */
    int  *n,             /* OUT - Vector length */
    MPI_Comm comm)       /* IN - Communicator */
{
    int datum_size;     /* Bytes per element */
    int i;
    FILE *infileptr;    /* Input file pointer */
    int local_els;      /* Elements on this proc */
    MPI_Status status;  /* Result of receive */
    int id;             /* Process rank */
    int p;              /* Number of processes */
    int x;              /* Result of read */ 
    datum_size = get_size (dtype); 
    MPI_Comm_size(comm, &p); 
    MPI_Comm_rank(comm, &id); 

    /* Process p-1 opens file, determines number of vector 
    elements, and broadcasts this value to the other processes. */ 
    if (id == (p-1)) {
        infileptr = fopen (s, "r"); 
        if (infileptr == NULL) *n = 0; 
        else fread (n, sizeof(int), 1, infileptr); 
    }
    MPI_Bcast (n, 1, MPI_INT, p-1, comm); 
    if (! *n) { 
        if (!id) { 
            printf ("Input file '%s' cannot be opened\n", s); 
            fflush (stdout); 
        }
    } 
    /* Block mapping of vector elements to processes */ 
    local_els = BLOCK_SIZE(id,p,*n); 

    /* Dynamically allocate vector. */ 
    *v = my_malloc (id, local_els * datum_size); 
    if (id == (p-1)) { 
        for (i = 0; i < p-1; i++) { 
            x = fread (*v, datum_size, BLOCK_SIZE(i,p,*n), infileptr); 
            MPI_Send (*v, BLOCK_SIZE(i,p,*n), dtype, i, DATA_MSG, comm); 
        }
        x = fread (*v, datum_size, BLOCK_SIZE(id,p,*n), infileptr); 
        fclose (infileptr); 
    } else { 
        MPI_Recv (*v, BLOCK_SIZE(id,p,*n), dtype, p-1, DATA_MSG, comm, &status); 
    } 
} 

/* Open a file containing a vector, read its contents, 
and replicate the vector among all processes in a communicator. */ 

void read_replicated_vector(
    char *s,             /* IN - File name */
    void **v,            /* OUT - Vector */
    MPI_Datatype dtype,  /* IN - Vector type */
    int  *n,             /* OUT - Vector length */
    MPI_Comm comm)       /* IN - Communicator */
{ 
    int datum_size;      /* Bytes per vector element */ 
    int i;
    int id;              /* Process rank */
    FILE *infileptr;     /* Input file pointer */
    int p;               /* Number of processes */ 

  
    MPI_Comm_rank (comm, &id); 
    MPI_Comm_size (comm, &p); 
    datum_size = get_size (dtype); 
    if (id == (p-1)) {
        infileptr = fopen (s, "r"); 
        if (infileptr == NULL) *n = 0; 
        else fscanf (infileptr, "%d ", n); 
        printf("n = %d\n", *n);
    } 
    MPI_Bcast (n, 1, MPI_INT, p-1, MPI_COMM_WORLD); 
    if (! *n) terminate (id, "Cannot open vector file"); 
    *v = my_malloc (id, *n * datum_size); 
    if (id == (p-1)) { 
        for (i = 0; i < *n; ++i)
            fscanf(infileptr, "%lf ", (double*)*v + i);
        fclose (infileptr); 
    }
    MPI_Bcast (*v, *n, dtype, p-1, MPI_COMM_WORLD); 
}

/******************** OUTPUT FUNCTIONS ********************/ 
/* 
* Print elements of a doubly subscripted array. 
*/ 
void print_submatrix ( 
    void            **a,        /* OUT - Doubly subscripted array */
    MPI_Datatype    dtype,      /* OUT - Type of array elements */
    int             rows,       /* OUT - Matrix rows */
    int             cols)       /* OUT - Matrix cols */ 
{ 
    int i, j; 
   
    for (i = 0; i < rows; i++) { 
        for (j = 0; j < cols; j++) { 

            if (dtype == MPI_DOUBLE) 
                printf ("%6.3f ", ((double **)a)[i][i]); 
            else {
                if (dtype == MPI_FLOAT) 
                printf ("%6.3f ", ((float **)a)[i][j]); 
                else if (dtype == MPI_INT) 
                printf ("%6d ", ((int **)a)[i][j]); 
            } 
        } 
        putchar ('\n'); 
    }
} 

/* • Print elements of a singly subscripted array. 
*/ 
void print_subvector (
    void            *a,         /* IN - Array pointer */ 
    MPI_Datatype    dtype,      /* IN - Array type */ 
    int             n)          /* IN - Array size */ 
{ 
    int i; 
    for (i = 0; i < n; i++) { 
        if (dtype == MPI_DOUBLE) 
            printf ("%6.3f ", ((double *)a)[i]); 
        else { 
            if (dtype == MPI_FLOAT) 
                printf ("%6.3f ", ((float *)a)[i]); 
            else if (dtype == MPI_INT) 
                printf ("%6d ", ((int *)a)[i]); 
        }
    }
} 

/* • Print a matrix distributed checkerboard fashion among • the processes in a communicator. 
*/ 

void print_checkerboard_matrix (
    void            **a,        /* IN -2D matrix */
    MPI_Datatype    dtype,      /* IN -Matrix element type */ 
    int             m,          /* IN -Matrix rows */ 
    int             n,          /* IN -Matrix columns */ 
    MPI_Comm        grid_comm)  /* IN - Communicator */ 
{
    void        *buffer;            /* Room to hold 1 matrix row */ 
    int         coords[2];          /* Grid coords of process sending elements */ 
    int         datum_size;         /* Bytes per matrix element  */
    int         els;                /* Elements received */ 
    int         grid_coord[2];     /* Coords of this process */ 
    int         grid_id;            /* Process rank in grid */ 
    int         grid_period[2];     /* Wraparound */ 
    int         grid_size[2];       /* Dims of process grid */
    int         i, j, k;             
    void        *laddr;             /* Where to put subrow */ 
    int         local_cols;         /* Matrix cols on this proc */ 
    int         p;                  /* Number of processes */ 
    int         src;                /* ID of proc with subrow */ 
    MPI_Status  status;             /* Result of receive */ 

    MPI_Comm_rank (grid_comm, &grid_id); 
    MPI_Comm_size (grid_comm, &p); 
    datum_size = get_size (dtype); 

    MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coord); 
    local_cols = BLOCK_SIZE(grid_coord[1], grid_size[1], n); 

    if (!grid_id) 
        buffer = my_malloc (grid_id, n * datum_size); 

    /* For each row of the process grid */ 
    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = i; 
        /* For each matrix row controlled by the process row */ 
        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], m); j++) {

            /* Collect the matrix row on grid process 0 and print it. */ 
            if (!grid_id) {
                for (k = 0; k < grid_size[1]; k++) {
                    coords[1] = k; 
                    MPI_Cart_rank (grid_comm, coords, &src); 
                    els = BLOCK_SIZE(k, grid_size[1], n); 
                    laddr = buffer + BLOCK_LOW(k, grid_size[1], n) * datum_size; 
                    if (src == 0) {
                        memcpy (laddr, a[j], els * datum_size); 
                    } else { 
                        MPI_Recv(laddr, els, dtype, src, 0, grid_comm, &status); 
                    } 
                }
                print_subvector (buffer, dtype, n); 
                putchar ('\n'); 
            } else if (grid_coord[0] == i) { 
                MPI_Send (a[j], local_cols, dtype, 0, 0, grid_comm); 
            }
        } 
    } 

    if (!grid_id) {
        free (buffer); 
        putchar ('\n'); 
    }
} 

/* 
Print a matrix that has a columnwise-block-striped data decomposition 
among the elements of a communicator. 
*/

void print_col_striped_matrix (
    void            **a,        /* IN - 2D array */ 
    MPI_Datatype    dtype,      /* IN - Type of matrix elements */ 
    int             m,          /* IN - Matrix rows */ 
    int             n,          /* IN - Matrix cols */ 
    MPI_Comm        comm)       /* IN - Communicator */ 
{
    MPI_Status      status;     /* Result of receive */ 
    int             datum_size; /* Bytes per matrix element */ 
    void            *buffer;    /* Enough room to hold 1 row */
    int             i, j;  
    int             id;         /* Process rank */ 
    int             p;          /* Number of processes */ 
    int*            rec_count;  /* Elements received per proc */ 
    int*            rec_disp;   /* Offset of each proc's block */ 

    MPI_Comm_rank (comm, &id); 
    MPI_Comm_size (comm, &p); 
    datum_size = get_size (dtype); 
    create_mixed_xfer_arrays (id, p, n, &rec_count, &rec_disp); 
    if (!id) 
        buffer = my_malloc (id, n*datum_size); 

    for (i = 0; i < m; i++) {
        MPI_Gatherv (a[i], BLOCK_SIZE(id,p,n), dtype, buffer, rec_count, 
                    rec_disp, dtype, 0, MPI_COMM_WORLD); 
        if (!id) { 
            print_subvector (buffer, dtype, n); 
            putchar ('\n'); 
        } 
    } 
    free (rec_count); 
    free (rec_disp); 
    if (!id) { 
        free (buffer); 
        putchar ('\n'); 
    } 
}

/* Print a matrix that is distributed in row-striped 
fashion among the processes in a communicator. */
void print_row_striped_matrix (
    void            **a,            /* IN - 2D array */ 
    MPI_Datatype    dtype,          /* IN - Matrix element type */ 
    int             m,              /* IN - Matrix rows */ 
    int             n,              /* IN - Matrix cols */ 
    MPI_Comm        comm)           /* IN - Communicator */ 
{ 
    MPI_Status      status;         /* Result of receive */
    void            *bstorage;      /* Elements received from another process */
    void            **b;            /* 2D array indexing into 'bstorage' */ 
    int             datum_size;     /* Bytes per element */ 
    int             i;              
    int             id;             /* Process rank */ 
    int             local_rows;     /* This proc's rows */ 
    int             max_block_size; /* Most matrix rows held by any process */ 
    int             prompt;         /* Dummy variable */ 
    int             p;              /* Number of processes */ 

    MPI_Comm_rank (comm, &id); 
    MPI_Comm_size (comm, &p); 
    local_rows = BLOCK_SIZE(id,p,m); 
    if (!id) { 
        print_submatrix (a, dtype, local_rows, n); 
        if (p > 1) { 
            datum_size = get_size (dtype); 
            max_block_size = BLOCK_SIZE(p-1,p,m); 
            bstorage = my_malloc (id, max_block_size * n * datum_size); 
            b = (void **) my_malloc (id, max_block_size * datum_size); 
            b[0] = bstorage; 
            for (i = 1; i < max_block_size; i++) { 
                b[i] = b[i-1] + n * datum_size; 
            } 
            for (i = 1; i < p; i++) { 
                MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG, MPI_COMM_WORLD); 
                MPI_Recv (bstorage, BLOCK_SIZE(i,p,m) * n, dtype, i,  RESPONSE_MSG, MPI_COMM_WORLD, &status); 
                print_submatrix (b, dtype, BLOCK_SIZE(i,p,m), n); 
            }
            free (b); 
            free (bstorage); 
        } 
        putchar ('\n'); 
    } else { 
        MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, MPI_COMM_WORLD, &status); 
        MPI_Send (*a, local_rows * n, dtype, 0, RESPONSE_MSG, MPI_COMM_WORLD); 
    } 
} 


/* 
Print a vector that is block distributed among the * processes in a communicator. 
*/
void print_block_vector (
    void            *v,     /* IN - Address of vector */ 
    MPI_Datatype    dtype,  /* IN - Vector element type */ 
    int             n,      /* IN - Elements in vector */ 
    MPI_Comm        comm)   /* IN - Communicator */ 
{
    int             datum_size;         /* Bytes per vector element */ 
    int             i; 
    int             prompt;             /* Dummy variable */ 
    MPI_Status      status;             /* Result of receive */ 
    void            *tmp;               /* Other process's subvector */ 
    int             id;                 /* Process rank */ 
    int             p;                  /* Number of processes */ 

    MPI_Comm_size (comm, &p); 
    MPI_Comm_rank (comm, &id); 
    datum_size = get_size (dtype); 
    if (!id) { 
        print_subvector (v, dtype, BLOCK_SIZE(id,p,n)); 
        if (p > 1) { 
            tmp = my_malloc (id, BLOCK_SIZE(p-1,p,n) * datum_size); 
            for (i = 1; i < p; i++) { 
                MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG, comm); 
                MPI_Recv (tmp, BLOCK_SIZE(i,p,n), dtype, i, RESPONSE_MSG, comm, &status); 
                print_subvector (tmp, dtype, BLOCK_SIZE(i,p,n)); 
            }
            free (tmp); 
        }
        printf ("\n\n"); 
     } else { 
         MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm, &status); 
         MPI_Send (v, BLOCK_SIZE(id,p,n), dtype, 0, RESPONSE_MSG, comm); 
     }
} 

/* 
* Print a vector that is replicated among the processes 
* in a communicator. 
*/

void print_replicated_vector (
    void            *v,             /* IN - Address of vector */ 
    MPI_Datatype    dtype,          /* IN - Vector element type */ 
    int             n,              /* IN - Elements in vector */ 
    MPI_Comm        comm)           /* IN - Communicator */ 
{ 
    int             id;             /* Process rank */ 

    MPI_Comm_rank (comm, &id); 

    if (!id) { 
        print_subvector (v, dtype, n); 
        printf ("\n\n");  
    } 
}  

void conv2d(
    D_TYPE **subs,
    D_TYPE *storage,
    int m,
    int n,
    D_TYPE **kernel_subs,
    D_TYPE *kernel_storage,
    int k,
    D_TYPE ***conv_subs,
    D_TYPE **conv_storage,
    MPI_Comm cart_comm
) {
    int rows, cols;
    D_TYPE **local_subs;
    D_TYPE *local_storage;
    int local_rows, local_cols; /* rows and cols for local_subs */
    int conv_rows, conv_cols;   /* rows and cols for conv_subs */
    int p;
    int grid_id;
    int grid_size[2];
    int grid_coords[2];
    int grid_periods[2];
    int coords[2];
    MPI_Status status;

    int dest_id, src_id;
    int i, j, i1, j1;

    MPI_Comm_size(cart_comm, &p);
    MPI_Comm_rank(cart_comm, &grid_id);
    MPI_Cart_get(cart_comm, 2, grid_size, grid_periods, grid_coords);

    rows = BLOCK_SIZE(grid_coords[0], grid_size[0], m);
    cols = BLOCK_SIZE(grid_coords[1], grid_size[1], n);

    local_rows = rows;
    local_cols = cols;
    if (grid_coords[1] != grid_size[1]-1)
        local_cols += (k - 1);
    if (grid_coords[0] != grid_size[0]-1)
        local_rows += (k - 1);
    
    // printf("[proc %d] rows %d cols %d local_rows %d local_cols %d\n", grid_id, rows, cols, local_rows, local_cols);

    local_storage = (D_TYPE *)my_malloc(grid_id, local_rows*local_cols*sizeof(D_TYPE));
    local_subs = (D_TYPE **)my_malloc(grid_id, local_rows*PTR_SIZE);
    for (i = 0; i < local_rows; ++ i)
        local_subs[i] = local_storage + i * local_cols;

    /* Initialization step 0: copy the rows x cols matrix to local_subs */
    for (i = 0; i < rows; ++ i)
        memcpy(local_subs[i], subs[i], cols*sizeof(D_TYPE));

    /* Initialization step 1: send the leftmost k-1 columns to the left */
    coords[0] = grid_coords[0];

    for (j = 0; j < rows; ++ j) {

        if (grid_coords[1] == 0) {
            if (grid_size[1] > 1) {
                coords[1] = grid_coords[1] + 1;
                MPI_Cart_rank(cart_comm, coords, &src_id);
                MPI_Recv(
                    local_subs[j]+cols, k-1, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
            }

        } else if (grid_coords[1] == grid_size[1]-1) {
            if (grid_size[1] > 1) {
                coords[1] = grid_coords[1] - 1;
                MPI_Cart_rank(cart_comm, coords, &dest_id);
                MPI_Send(
                    subs[j], k-1, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);    
            }

        } else {
            coords[1] = grid_coords[1] - 1;
            MPI_Cart_rank(cart_comm, coords, &dest_id);
            coords[1] = grid_coords[1] + 1;
            MPI_Cart_rank(cart_comm, coords, &src_id);

            if (grid_coords[1] % 2) {
                MPI_Send(
                    subs[j], k-1, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);
                MPI_Recv(
                    local_subs[j]+cols, k-1, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
            } else {
                MPI_Recv(
                    local_subs[j]+cols, k-1, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
                MPI_Send(
                    subs[j], k-1, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);
            }
        }
    }

    /* Initialization step 2: send the upmost k-1 rows to the up */
    coords[1] = grid_coords[1];

    for (j = 0; j < k-1; ++ j) {

        if (grid_coords[0] == 0) {
            if (grid_size[0] > 1) {
                coords[0] = grid_coords[0] + 1;
                MPI_Cart_rank(cart_comm, coords, &src_id);
                MPI_Recv(
                    local_subs[rows+j], cols, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
            }

        } else if (grid_coords[0] == grid_size[0]-1) {
            if (grid_size[0] > 1) {
                coords[0] = grid_coords[0] - 1;
                MPI_Cart_rank(cart_comm, coords, &dest_id);
                MPI_Send(
                    subs[j], cols, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);    
            }

        } else {
            coords[0] = grid_coords[0] - 1;
            MPI_Cart_rank(cart_comm, coords, &dest_id);
            coords[0] = grid_coords[0] + 1;
            MPI_Cart_rank(cart_comm, coords, &src_id);

            if (grid_coords[0] % 2) {
                MPI_Send(
                    subs[j], cols, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);
                MPI_Recv(
                    local_subs[rows+j], cols, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
            } else {
                MPI_Recv(
                    local_subs[rows+j], cols, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
                MPI_Send(
                    subs[j], cols, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);
            }
        }
    }

    /* Initialization step 3: send the up-left (k-1)x(k-1) sub-matrix to the up-left */
    for (j = 0; j < k-1; ++ j) {

        if ((grid_coords[0] == 0 && grid_coords[1] == grid_size[1]-1) || 
                (grid_coords[0] == grid_size[0]-1 && grid_coords[1] == 0)) {
            ;
        } else if (grid_coords[0] == 0 || grid_coords[1] == 0) {
            if (grid_size[0] > 1 && grid_size[1] > 1) {
                coords[0] = grid_coords[0] + 1;
                coords[1] = grid_coords[1] + 1;
                MPI_Cart_rank(cart_comm, coords, &src_id);
                MPI_Recv(
                    &local_subs[rows+j][cols], k-1, MPI_TYPE, src_id, 
                    CONV_INIT_TAG, cart_comm, &status);
            }

        } else if (grid_coords[0] == grid_size[0]-1 || grid_coords[1] == grid_size[1]-1) {
            if (grid_size[0] > 1 && grid_size[1] > 1) {
                coords[0] = grid_coords[0] - 1;
                coords[1] = grid_coords[1] - 1;
                MPI_Cart_rank(cart_comm, coords, &dest_id);
                MPI_Send(
                    subs[j], k-1, MPI_TYPE, dest_id, 
                    CONV_INIT_TAG, cart_comm);
            }
            
        } else {
            coords[0] = grid_coords[0] + 1;
            coords[1] = grid_coords[1] + 1;
            MPI_Cart_rank(cart_comm, coords, &src_id);
            coords[0] = grid_coords[0] - 1;
            coords[1] = grid_coords[1] - 1;
            MPI_Cart_rank(cart_comm, coords, &dest_id);

            // Send on receive
            MPI_Recv(
                &local_subs[rows+j][cols], k-1, MPI_TYPE, src_id, 
                CONV_INIT_TAG, cart_comm, &status);
            MPI_Send(
                subs[j], k-1, MPI_TYPE, dest_id, 
                CONV_INIT_TAG, cart_comm);
        }
    }
    if (local_rows >= k && local_cols >= k) {
        conv_rows = local_rows - k + 1;
        conv_cols = local_cols - k + 1;
        *conv_storage = (D_TYPE *)my_malloc(grid_id, conv_rows*conv_cols*sizeof(D_TYPE));
        *conv_subs = (D_TYPE **)my_malloc(grid_id, conv_rows*PTR_SIZE);
        for (i = 0; i < conv_rows; ++ i) {
            (*conv_subs)[i] = *conv_storage + i * conv_cols;
        }

        for (i = 0; i < conv_rows; ++ i) {
            for (j = 0; j < conv_cols; ++ j) {
                 (*conv_subs)[i][j] = 0.0;
                for (i1 = 0; i1 < k; ++ i1) {
                    for (j1 = 0; j1 < k; ++ j1) {
                        (*conv_subs)[i][j] += local_subs[i+i1][j+j1] * kernel_subs[i1][j1];
                    }
                }
            }
        }

    } else {
        *conv_storage = NULL;
        *conv_subs = NULL;
    }

    free(local_storage);
    free(local_subs);
}


void print_conv_result(
    D_TYPE **conv_subs,
    int m,
    int n,
    int kernel_size,
    MPI_Datatype dtype,
    MPI_Comm cart_comm
) {
    int grid_id;
    int p;
    int grid_size[2];
    int grid_periods[2];
    int grid_coords[2];
    int local_rows, local_cols;     /* for sub-matrix */
    int conv_rows, conv_cols;       /* for conv_subs */
    int coords[2];
    int src_id;
    int i, j, t;
    int k = kernel_size;
    D_TYPE *row_buffer;
    int row_buffer_idx;
    int *conv_col_array;
    MPI_Status status;

    MPI_Comm_size(cart_comm, &p);
    MPI_Comm_rank(cart_comm, &grid_id);
    MPI_Cart_get(cart_comm, 2, grid_size, grid_periods, grid_coords);
    
    if (grid_id == 0) {
        row_buffer = (D_TYPE *)my_malloc(grid_id, (n-k+1)*sizeof(D_TYPE));
    }

    conv_col_array = (int *)my_malloc(grid_id, grid_size[1]*sizeof(int));
    for (i = 0; i < grid_size[1]; ++ i) {
        local_cols = BLOCK_SIZE(i, grid_size[1], n);
        if (i != grid_size[1]-1)
            local_cols += (k - 1);
        conv_col_array[i] = MAX(local_cols-k+1, 0);
    }

    for (i = 0; i < grid_size[0]; ++ i) {
        coords[0] = i;

        local_rows = BLOCK_SIZE(i, grid_size[0], m);
        if (i != grid_size[0]-1)
            local_rows += (k - 1);
        conv_rows = MAX(local_rows-k+1, 0);

        if (conv_rows == 0) continue;

        for (j = 0; j < conv_rows; ++ j) {
            row_buffer_idx = 0;

            for (t = 0; t < grid_size[1]; ++ t) {                    
                conv_cols = conv_col_array[t];
                if (conv_cols == 0) continue;
                
                coords[1] = t;
                MPI_Cart_rank(cart_comm, coords, &src_id);

                if (src_id == 0) {
                    if (grid_id == 0) {
                        memcpy(
                            row_buffer+row_buffer_idx, conv_subs[j], 
                            conv_cols*sizeof(D_TYPE));
                    }
                } else {
                    if (grid_id == 0) {
                        MPI_Recv(
                            row_buffer+row_buffer_idx, conv_cols, dtype, src_id, 
                            CONV_PRINT_TAG, cart_comm, &status);
                    } else if (grid_id == src_id) {
                        MPI_Send(
                            conv_subs[j], conv_cols, dtype, 0, 
                            CONV_PRINT_TAG, cart_comm);
                    }
                }
                row_buffer_idx += conv_cols;
            }

            if (grid_id == 0) {
                for (t = 0; t < n-k+1; ++ t) {
                    printf(PRINT_FORMAT, row_buffer[t]);
                }
                printf("\n");
            }
        }
    }

    free(conv_col_array);
    if (grid_id == 0)
        free(row_buffer);
}

int matmul(
    D_TYPE **A_subs, D_TYPE *A_storage, int m1, int n1,
    D_TYPE **B_subs, D_TYPE *B_storage, int m2, int n2,
    D_TYPE ***C_subs, D_TYPE **C_storage,
    MPI_Datatype dtype,
    MPI_Comm cart_comm
) {
    int id, p, size[2], periods[2], coords[2], local_coords[2];
    int src_id, dest_id, datum_size;
    D_TYPE *recv_buffer;
    int A_rows, A_cols, B_rows, B_cols, C_rows, C_cols, i, j, k, a;
    MPI_Status status;
    MPI_Comm_size(cart_comm, &p);
    MPI_Comm_rank(cart_comm, &id);
    MPI_Cart_get(cart_comm, 2, size, periods, local_coords);
    MPI_Type_size(dtype, &datum_size);
    if (n1 != m2) {
        printf("[ERROR id=%d] Matrix A (%d, %d) Matrix B (%d, %d)\n", 
            id, m1, n1, m2, n2);
        MPI_Abort(MPI_COMM_WORLD, -6);} 
    A_rows = BLOCK_SIZE(local_coords[0], size[0], m1);
    A_cols = BLOCK_SIZE(local_coords[1], size[1], n1);
    B_rows = BLOCK_SIZE(local_coords[0], size[0], m2);
    B_cols = BLOCK_SIZE(local_coords[1], size[1], n2);
    C_rows = A_rows;
    C_cols = B_cols;
    (*C_storage) = (D_TYPE *)my_malloc(id, datum_size*C_rows*C_cols);
    (*C_subs) = (D_TYPE **)my_malloc(id, sizeof(D_TYPE*)*C_rows);
    for (i = 0; i < C_rows; ++ i) {(*C_subs)[i] = (*C_storage) + i * C_cols;}
    for (i = 0; i < C_rows; ++ i) {for (j = 0; j < C_cols; ++ j) {(*C_subs)[i][j] = 0.0;}}
    recv_buffer = (D_TYPE *)my_malloc(
        id, MAX(A_rows,B_rows)*MAX(A_cols,B_cols)*datum_size);
    if (local_coords[0]) {
        coords[0] = local_coords[0];
        coords[1] = MOD(local_coords[1] - local_coords[0], size[1]);
        MPI_Cart_rank(cart_comm, coords, &dest_id);
        coords[1] = MOD(local_coords[1] + local_coords[0], size[1]);
        MPI_Cart_rank(cart_comm, coords, &src_id);
        if (local_coords[1] == 0) {
            MPI_Send(A_storage, A_rows*A_cols, dtype, dest_id, INITIAL_TAG, cart_comm);
            MPI_Recv(recv_buffer, A_rows*A_cols, dtype, src_id,  INITIAL_TAG, cart_comm, &status);
        } else {MPI_Recv(recv_buffer, A_rows*A_cols, dtype, src_id, INITIAL_TAG, cart_comm, &status);
            MPI_Send( A_storage, A_rows*A_cols, dtype, dest_id, INITIAL_TAG, cart_comm);
        }
        memcpy(A_storage, recv_buffer, A_rows*A_cols*sizeof(D_TYPE));
    }
    if (local_coords[1]) {
        coords[1] = local_coords[1];
        coords[0] = MOD(local_coords[0] - local_coords[1], size[0]);
        MPI_Cart_rank(cart_comm, coords, &dest_id);
        coords[0] = MOD(local_coords[0] + local_coords[1], size[0]);
        MPI_Cart_rank(cart_comm, coords, &src_id);
        if (local_coords[0] == 0) {
            MPI_Send(B_storage, B_rows*B_cols, dtype, dest_id, INITIAL_TAG, cart_comm);
            MPI_Recv(recv_buffer, B_rows*B_cols, dtype, src_id, INITIAL_TAG, cart_comm, &status);
        } else {MPI_Recv(recv_buffer, B_rows*B_cols, dtype, src_id, INITIAL_TAG, cart_comm, &status);
            MPI_Send(B_storage, B_rows*B_cols, dtype, dest_id, INITIAL_TAG, cart_comm);}
        memcpy(B_storage, recv_buffer, B_rows*B_cols*datum_size);}
    if (PRINT_FLAG) {
        if (id == 0) {
            printf(">> After alignment (initialization)\n");
            printf(">> Aligned factor matrix A\n");
        }
        print_checkerboard_matrix((void **)A_subs, dtype, m1, n1, cart_comm);
        if (id == 0) printf("\n");

        if (id == 0) printf(">> Aligned factor matrix B\n");
        print_checkerboard_matrix((void **)B_subs, dtype, m2, n2, cart_comm);
        if (id == 0) printf("\n");
    }
    for (k = 0; k < size[0]; ++ k) {
        for (i = 0; i < C_rows; ++ i) {
            for (j = 0; j < C_cols; ++ j) {
                for (a = 0; a < A_cols; ++ a) {
                    (*C_subs)[i][j] += (A_subs[i][a] * B_subs[a][j]);}}}
        if (k != size[0] - 1) {
            coords[0] = local_coords[0];
            coords[1] = MOD(local_coords[1] - 1, size[1]);
            MPI_Cart_rank(cart_comm, coords, &dest_id);
            coords[1] = MOD(local_coords[1] + 1, size[1]);
            MPI_Cart_rank(cart_comm, coords, &src_id);
            if (local_coords[1] % 2) {
                MPI_Send(A_storage, A_rows*A_cols, dtype, dest_id, CANNON_TAG, cart_comm);
                MPI_Recv(recv_buffer, A_rows*A_cols, dtype, src_id, CANNON_TAG, cart_comm, &status);
            } else {
                MPI_Recv(recv_buffer, A_rows*A_cols, dtype, src_id, CANNON_TAG, cart_comm, &status);
                MPI_Send(A_storage, A_rows*A_cols, dtype, dest_id, CANNON_TAG, cart_comm);
            }
            memcpy(A_storage, recv_buffer, A_rows*A_cols*datum_size);
            // coords[1] = local_coords[1];
            // coords[0] = MOD(local_coords[0] - 1, size[0]);
            // MPI_Cart_rank(cart_comm, coords, &dest_id);
            coords[1] = local_coords[1];
            coords[0] = MOD(local_coords[0] - 1, size[0]);
            MPI_Cart_rank(cart_comm, coords, &dest_id);
            coords[0] = MOD(local_coords[0] + 1, size[0]);
            MPI_Cart_rank(cart_comm, coords, &src_id);
            if (local_coords[0] % 2) {
                MPI_Send(B_storage, B_rows*B_cols, dtype, dest_id, CANNON_TAG, cart_comm);
                MPI_Recv(recv_buffer, B_rows*B_cols, dtype, src_id, CANNON_TAG, cart_comm, &status);
            } else {
                MPI_Recv(recv_buffer, B_rows*B_cols, dtype, src_id, CANNON_TAG, cart_comm, &status);
                MPI_Send(B_storage, B_rows*B_cols, dtype, dest_id, CANNON_TAG, cart_comm);
            }
            memcpy(B_storage, recv_buffer, B_rows*B_cols*datum_size);}}
    free(recv_buffer);
}
