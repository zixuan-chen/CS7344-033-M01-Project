#include "MyMPI.h"

int main(int argc, char *argv[]) {
    int p, id;
    D_TYPE **subs,*storage,**kernel_subs,*kernel_storage,**conv_subs,*conv_storage;
    int m, n, k, k2, i, j;
    MPI_Comm cart_comm;
    int size[2],periods[2],local_coords[2];
    double conv_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (argc < 3) {if (id == 0)
        printf("ERROR Missing arguments. Note:\
            %s <matrix path> <kernel path>\n", argv[0]);
        MPI_Finalize();
        return 0;}
    size[0] = size[1] = 0;
    MPI_Dims_create(p, 2, size);
    if (size[0] != size[1]) {
        printf("ERROR dimension must be squired(e.g. 4, 16, 25)\n");
        MPI_Abort(MPI_COMM_WORLD, NODE_NOT_SQAURE_ERROR);}
    periods[0] = periods[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, id, 2, local_coords);
    read_checkerboard_matrix(
        argv[1], (void ***)&subs, (void **)&storage, 
        MPI_TYPE, &m, &n, cart_comm);
    if (PRINT_FLAG) {
        if (id == 0) printf(">> Feature matrix\n");
        print_checkerboard_matrix((void **)subs, MPI_TYPE, m, n, cart_comm);
        if (id == 0) printf("\n");}
    read_simple_matrix(
        argv[2], (void ***)&kernel_subs, (void **)&kernel_storage, 
        MPI_TYPE, &k, &k2, MPI_COMM_WORLD);
    if (k != k2) {
        printf("[error] Kernel must take the square form\n");
        MPI_Abort(MPI_COMM_WORLD, MATRIX_SIZE_MISMATCH_ERROR);
    }
    if (PRINT_FLAG) {
        if (id == 0) {
            printf(">> Kernel matrix\n");
            for (i = 0; i < k; ++ i) {
                for (j = 0; j < k; ++ j) {
                    printf(PRINT_FORMAT, kernel_subs[i][j]);}
                printf("\n");}
            printf("\n");}}
    MPI_Barrier(MPI_COMM_WORLD);
    conv_time = -MPI_Wtime();
    conv2d(subs, storage, m, n, kernel_subs,
    kernel_storage, k, &conv_subs, &conv_storage, cart_comm);
    MPI_Barrier(MPI_COMM_WORLD);
    conv_time += MPI_Wtime();

    if (PRINT_FLAG) {
        if (id == 0) printf("[result] Convolution operation result\n");
        print_conv_result(conv_subs, m, n, k, MPI_TYPE, cart_comm);
        if (id == 0) printf("\n");
    }

    if (id == 0) 
        printf(
            "Matrix size: (%d, %d), Kernel size (%d, %d). Total elapsed time: %10.6f\n", 
            m, n, k, k2, conv_time);

    MPI_Finalize();

    return 0;
}