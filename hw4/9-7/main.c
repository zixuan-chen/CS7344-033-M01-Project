#include "mpi.h"
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <sys/stat.h> 
#include <ftw.h> 

/******** ********* ******** MACROS ************* ***** ********/ 
#define DATA_MSG            0 
#define PROMPT_MSG          1 
#define RESPONSE_MSG        2 
#define OPEN_FILE_ERROR    -1 
#define MALLOC_ERROR       -2 
#define TYPE_ERROR         -3 
#define MIN(a,b)            ((a)<(b)?(a):(b)) 
#define BLOCK_LOW(id,p,n)   ((id) * (n)/(p)) 
#define BLOCK_HIGH(id,p,n)  (BLOCK_LOW((id)+1,p,n)-1) 
#define BLOCK_SIZE(id,p,n)  (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1) 
#define BLOCK_OWNER(j,p,n)  (((p) * ((j)+1)-1)/(n)) 
#define PTR_SIZE            (sizeof(void*)) 
#define CEILING(i,j)        (((i)+(j)-1)/(j)) 
#define DICT_SIZE_MSG       0 /* Msg has dictionary size */ 
#define FILE_NAME_MSG       1 /* Msg is file name */ 
#define VECTOR_MSG          2 /* Msg is profile */ 
#define EMPTY_MSG           3 /* Msg is empty */ 
#define INPUT_ARG           1 /* matrix argument */ 
#define OUTPUT_ARG          2 /* vector argument */ 

typedef unsigned char uchar; 

int main (int argc, char *argv[]) { 
    int id; /* Process rank */ 
    int p; /* Number of processes */ 
    MPI_Comm worker_comm; /* Workers-only communicator  */

    void manager (int, char **, int); 
    void worker (int, char **, MPI_Comm); 

    MPI_Init (&argc, &argv); 
    MPI_Comm_rank (MPI_COMM_WORLD, &id); 
    MPI_Comm_size (MPI_COMM_WORLD, &p); 
    if (argc != 2) { 
        if (!id) {
            printf ("Program needs one arguments:\n"); 
            printf ("%s <in_path>\n", argv[0]);
        }
    } else if (p < 1) { 
        printf ("Program needs at least two processes\n"); 
    } else {
        if (!id) { 
            MPI_Comm_split (MPI_COMM_WORLD, MPI_UNDEFINED, 
                            id, &worker_comm); 
            manager (argc, argv, p); 
        }else { 
            MPI_Comm_split (MPI_COMM_WORLD, 0, id, &worker_comm); 
            worker (argc, argv, worker_comm); 
        }
    }
    MPI_Finalize(); 
    return 0; 
}

void manager (int argc, char *argv[], int p) {

    double *a;         /* Store matrix here */
    double *b;          /* Store vector here */ 
    double *c;
    double res;
    int n;              /* column of matrix */ 
    int m;              /* row of matrix */
    int i, j; 
    int assign_cnt;
    // MPI_Request pending;/* Handle for recv request */ 
    int src;            /* Message source process */ 
    MPI_Status status;  /* Message status */ 
    int tag;            /* Message tag */ 
    int terminated;     /* Count of terminated procs */ 
    FILE * infileptr;
    /* read matrix and vectors */
    
    infileptr = fopen(argv[INPUT_ARG], "r");
    if (infileptr == NULL){
        m = 0; n = 0;
        printf("failed to open file: %s\n", argv[INPUT_ARG]);
        MPI_Abort(MPI_COMM_WORLD, OPEN_FILE_ERROR);
    }else{
        fscanf(infileptr, "%d %d", &m, &n);
        printf("m = %d, n = %d\n", m, n);
    }
    a = (double*) malloc(m * n * sizeof(double));
    b = (double*) malloc(n * sizeof(double));
    c = (double*) malloc(m *sizeof(double));
    if (a == NULL || b == NULL || c == NULL){
        printf("error in manager process: can allocate enough \
                memory for a, b, c\n");
        MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
    }
    for (i = 0; i < m; ++i){
        for (j = 0; j < n; ++j)
            fscanf(infileptr, "%lf", a + i * n + j);
    }
    for (j = 0; j < n; ++j)
        fscanf(infileptr, "%lf", b + j);
    fclose(infileptr);

    /* distribute copy of vector to all workers */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Respond to requests by workers. */
    terminated = 0; 
    assign_cnt = 0; 
    do {
        MPI_Recv (&res, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
                    MPI_COMM_WORLD, &status);
        src = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        if (tag == VECTOR_MSG){
            c[assign_cnt] = res;
            assign_cnt++;
        }
        /* Assign more work or tell worker to stop. */ 
        if (assign_cnt < m) { 
            MPI_Send (a + assign_cnt * n, n, MPI_DOUBLE, src, 
                        FILE_NAME_MSG, MPI_COMM_WORLD); 
        } else { 
            MPI_Send (NULL, 0, MPI_CHAR, src, FILE_NAME_MSG, 
                    MPI_COMM_WORLD); 
            terminated++; 
        } 
    } while (terminated < (p-1)); 
    
    printf("result: ");
    for (int i = 0; i < m; ++i){
        printf("%lf ", c[i]);
    }
    printf("\n");
    free(a);free(b);free(c);
} 

void worker (int argc, char *argv[], MPI_Comm worker_comm) { 
    double *b;
    double *row_a; 
    double res;
    int n; /* Profile vector size */ 
    int i; 
    int name_len; /* Chars in file name */ 
    MPI_Request pending; /* Handle for MPI_Isend */ 
    MPI_Status status; /* Info about message */ 
    int worker_id; /* Rank in worker_comm */ 

    MPI_Comm_rank(worker_comm, &worker_id);
    
    b = (double*) malloc(n * sizeof(double));
    row_a = (double*) malloc(n *sizeof(double));
    if(b == NULL || row_a == NULL)
    {
        if(!worker_id)
            printf("error in worker process: can't allocate \
            enough memory for b and row_a\n");
        MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    MPI_Isend(NULL, 0, MPI_DOUBLE, 0, EMPTY_MSG, MPI_COMM_WORLD, 
                &pending);
    for(;;){
        MPI_Probe(0, FILE_NAME_MSG, MPI_COMM_WORLD, &status);
        MPI_Get_count (&status, MPI_DOUBLE, &name_len);
        if (name_len != n) break;
        
        MPI_Recv(row_a, n, MPI_DOUBLE, 0, FILE_NAME_MSG,
                    MPI_COMM_WORLD, &status);

        res = 0;
        for (i = 0; i < n; ++i) res += row_a[i] * b[i];
        
        MPI_Send(&res, 1, MPI_DOUBLE, 0, VECTOR_MSG, MPI_COMM_WORLD);
    }
    free(b); free(row_a);
}


