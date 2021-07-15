#include "MyMPI.h"

int main(int argc, char *argv[]) {
    int p, id, size[2], periodic[2], local_coords[2], m1, n1, m2, n2;
    MPI_Comm cart_comm;
    double elapsed_time;
    D_TYPE **A_subs, *A_storage, **B_subs, *B_storage, **C_subs, *C_storage;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    size[0] = size[1] = 0;
    MPI_Dims_create(p, 2, size);
    if (size[0] != size[1]) {
        printf("[error] Only support squared nodes \n");
        MPI_Abort(MPI_COMM_WORLD, NODE_NOT_SQAURE_ERROR);}
    periodic[0] = periodic[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 1, &cart_comm);
    MPI_Cart_coords(cart_comm, id, 2, local_coords);
    read_checkerboard_matrix("matrix/A_1024.txt", (void ***)&A_subs, (void **)&A_storage, 
        MPI_TYPE, &m1, &n1, cart_comm);
    read_checkerboard_matrix("matrix/B_1024.txt", (void ***)&B_subs, (void **)&B_storage, 
        MPI_TYPE, &m2, &n2, cart_comm);
    if (PRINT_FLAG) {if (id == 0) printf(">>print matrix A\n");
        print_checkerboard_matrix((void **)A_subs, MPI_TYPE, m1, n1, cart_comm);
        if (id == 0) printf("\n");
        if (id == 0) printf(">>print matrix B\n");
        print_checkerboard_matrix((void **)B_subs, MPI_TYPE, m2, n2, cart_comm);
        if (id == 0)  printf("\n");}
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    matmul(A_subs, A_storage, m1, n1, B_subs, B_storage, m2, n2, 
        &C_subs, &C_storage, MPI_TYPE, cart_comm);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    if (PRINT_FLAG) {if (id == 0) printf("[result] Product matrix C\n");
        print_checkerboard_matrix((void **)C_subs, MPI_TYPE, m1, n2, cart_comm);
        if (id == 0) printf("\n");}
    if (id == 0) printf(
            "Multiplication size: (%d, %d) x (%d, %d). Total elapsed time: %10.6f\n", 
            m1, n1, m2, n2, elapsed_time);
    MPI_Finalize();
    return 0;
}