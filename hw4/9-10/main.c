#include "mpi.h"
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <sys/stat.h> 
#include <math.h>
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
#define SUCCESS_MSG         1 /* success in finding a perfect number */ 
#define FAIL_MSG            2 /* fail in finding a perfect number */ 
#define EMPTY_MSG           3 /* Msg is empty */ 
#define INPUT_ARG           1 /* matrix argument */ 
#define OUTPUT_ARG          2 /* vector argument */ 
#define NUM_PERFECT         8
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

    if (p < 1) { 
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

    int *e;         /* Store matrix here */
    int n;
    int i;
    int assign_cnt;
    int src;            /* Message source process */ 
    MPI_Status status;  /* Message status */ 
    int tag;            /* Message tag */ 
    int terminated;     /* Count of terminated procs */ 
    double perfect_number;
    
    e = (int *) malloc (NUM_PERFECT * sizeof(int));

    if (e == NULL){
        printf("error in manager process: can allocate enough \
                memory for e\n");
        MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
    }

    /* Respond to requests by workers. */
    terminated = 0; 
    assign_cnt = 0; 
    i = 1;
    do {
        MPI_Recv (&n, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
                    MPI_COMM_WORLD, &status);
        src = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        
        if (tag == SUCCESS_MSG){
            e[assign_cnt] = n;
            assign_cnt++;
        }
        /* Assign more work or tell worker to stop. */ 
        if (assign_cnt < NUM_PERFECT) { 
            MPI_Send (&i, 1, MPI_INT, src, SUCCESS_MSG, MPI_COMM_WORLD); 
        } else { 
            MPI_Send (NULL, 0, MPI_INT, src, SUCCESS_MSG, MPI_COMM_WORLD); 
            terminated++; 
        } 
        ++i;
    } while (terminated < (p-1)); 
    
    printf("the first %d perfect numbers : \n", NUM_PERFECT);
    for (int i = 0; i < NUM_PERFECT; ++i){
        perfect_number = (pow(2, e[i])-1) * pow(2, e[i]-1);
        printf("%d: %.0lf\n", i, perfect_number);
    }
    free(e);
} 

void worker (int argc, char *argv[], MPI_Comm worker_comm) { 
    int n; /* Profile vector size */ 
    int name_len; /* Chars in file name */ 
    MPI_Request pending; /* Handle for MPI_Isend */ 
    MPI_Status status; /* Info about message */ 
    int worker_id; /* Rank in worker_comm */ 

    int is_prime(int x);

    MPI_Comm_rank(worker_comm, &worker_id);

    MPI_Isend(NULL, 0, MPI_INT, 0, EMPTY_MSG, MPI_COMM_WORLD, &pending);
    for(;;){
        MPI_Probe(0, SUCCESS_MSG, MPI_COMM_WORLD, &status);
        MPI_Get_count (&status, MPI_INT, &name_len);
        if (!name_len) break;
        
        MPI_Recv(&n, 1, MPI_INT, 0, SUCCESS_MSG, MPI_COMM_WORLD, &status);

        if (is_prime(n)){
            MPI_Send(&n, 1, MPI_INT, 0, SUCCESS_MSG, MPI_COMM_WORLD);
        }else{
            MPI_Send(&n, 1, MPI_INT, 0, FAIL_MSG, MPI_COMM_WORLD);
        }
    }
}

int is_prime(int n){
    unsigned long long x = pow(2, n)-1;
    
    if((x < 2) || ((x%2 == 0) && (x!=2)))
        return 0;
    
    for(int i = 3; i < sqrt(x); i+=2){
        if (x%i == 0)
            return 0;
    }
    return 1;
}