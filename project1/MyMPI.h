/* MyMPI.h
*
* Header file for a library of matrix/vector 
* input/output/redistribution functions. 
• Programmed by Michael J. Quinn
• Last modification: 4 September 2002 
*/

/******** ********* ******** MACROS ************* ***** ********/ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define D_TYPE              double
#define MPI_TYPE            MPI_DOUBLE
#define PRINT_FORMAT        "%.1lf\t"
#define SCANF_FORMAT        "%lf"
#define DATA_MSG            0 
#define PROMPT_MSG          1 
#define RESPONSE_MSG        2 
#define OPEN_FILE_ERROR    -1 
#define MALLOC_ERROR       -2 
#define TYPE_ERROR         -3 
#define INITIAL_TAG         1
#define CANNON_TAG          2
#define PRINT_TAG           3
#define ARGS_ERROR          -5
#define MATRIX_SIZE_MISMATCH_ERROR -6
#define NODE_NOT_SQAURE_ERROR -7
#define CONV_INIT_TAG       1
#define CONV_PRINT_TAG      2
#ifndef PRINT_FLAG
#define PRINT_FLAG          0
#endif

#define MIN(a,b)            ((a)<(b)?(a):(b)) 
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MOD(x,m) (((x)+(m))%(m))
#define BLOCK_LOW(id,p,n)   ((id) * (n)/(p)) 
#define BLOCK_HIGH(id,p,n)  (BLOCK_LOW((id)+1,p,n)-1) 
#define BLOCK_SIZE(id,p,n)  (BLOCK_HIGH(id+1,p,n)-BLOCK_LOW(id,p,n)) 
#define BLOCK_OWNER(j,p,n)  (((p) * ((j)+1)-1)/(n)) 
#define PTR_SIZE            (sizeof(void*)) 
#define CEILING(i,j)        (((i)+(j)-1)/(j)) 


/********** ***** ** MISCELLANEOUS FUNCTIONS ******* ***** *****/ 
void terminate (int, char *); 
void* my_malloc (int, int);
/********** ***** DATA DISTRIBUTION FUNCTIONS ** ****** *******/ 
void replicate_block_vector (void *, int, void *, MPI_Datatype, MPI_Comm); 
void create_mixed_xfer_arrays (int, int, int, int**, int**); 
void create_uniform_xfer_arrays (int, int, int, int**,int**); 
void conv2d(D_TYPE **, D_TYPE *, int, int, D_TYPE **,D_TYPE *,int,D_TYPE ***,D_TYPE **,MPI_Comm);
/****************** INPUT FUNCTIONS ************* ******* ****/ 
void read_checkerboard_matrix (char *, void ***, void **, MPI_Datatype, int *, int *, MPI_Comm); 
void read_col_striped_matrix (char *, void ***, void **, MPI_Datatype, int *, int *, MPI_Comm); 
void read_row_striped_matrix (char *, void ***, void **, MPI_Datatype, int *, int *, MPI_Comm);
void read_block_vector (char *, void **, MPI_Datatype, int *, MPI_Comm); 
void read_replicated_vector (char *, void **, MPI_Datatype, int *, MPI_Comm); 
void read_simple_matrix(char *, void ***, void **, MPI_Datatype, int *, int *, MPI_Comm);

/*** ***** ***** ***** OUTPUT FUNCTIONS ******** ******** *******/ 
void print_checkerboard_matrix (void **, MPI_Datatype, int, int, MPI_Comm); 
void print_col_striped_matrix (void **, MPI_Datatype, int, int, MPI_Comm); 
void print_row_striped_matrix (void **, MPI_Datatype, int, int, MPI_Comm); 
void print_block_vector (void *, MPI_Datatype, int, MPI_Comm); 
void print_replicated_vector (void *, MPI_Datatype, int, MPI_Comm); 
void print_conv_result(D_TYPE **,int,int,int,MPI_Datatype,MPI_Comm);
int matmul(
    D_TYPE **, D_TYPE *, int, int,
    D_TYPE **, D_TYPE *, int , int ,
    D_TYPE ***, D_TYPE **,
    MPI_Datatype ,
    MPI_Comm);