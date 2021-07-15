SEP="----------------------------------------------------------"
NP=4
TEST_PROC=4 

echo Student name: Zixuan Chen
echo Student ID: 120033910005
echo Parallel Programming -- Project 1
echo


# Exercise 1: MPI_ALLGATHER
echo MPI_ALLGATHER
echo ${SEP}
mpirun -np ${NP} ./my_MPI_Allgather ${NP}
echo ${SEP}
echo


# Exercise 2: Gemm
echo Large matrix multiplication
echo ${SEP}
echo 'random generate 1024x1024 matrix...'
mpirun -np ${NP} ./matmul ./matrix/A_1024.txt ./matrix/B_1024.txt
echo ${SEP}
echo

# Exercise 2.2: convolution
echo Convolution operation
echo 'Using a 4 * 4 kernel to do the pooling operation for the matrix.'
mpirun -np ${NP} ./conv ./matrix/conv_1024.txt ./matrix/kernel.txt
echo ${SEP}
echo

# Exercise 2.3: pooling
echo Pooling operation
echo 'Using a 4 * 4 kernel to do the convolution operation for the matrix'
mpirun -np ${NP} ./conv ./matrix/pool_1024.txt ./matrix/pool_kernel_4.txt
echo ${SEP}
echo

# Exercise 3: Wordcount
echo Wordcount
echo ${SEP}
echo 'wordcount small files'
mpirun -np ${NP} ./wordcount_small
echo
echo 'wordcount big file'
mpirun -np ${NP} ./wordcount_big
echo ${SEP}
echo