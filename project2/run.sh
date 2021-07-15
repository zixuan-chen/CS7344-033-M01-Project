SEP="----------------------------------------------------------"

echo Student name: Zixuan Chen
echo Student ID: 120033910005
echo Parallel Programming -- Project 2
echo

echo SEP
echo Compiling Program
make clean
make
echo SEP

# Exercise 1: Monte Carlo
echo Exercise 1: Monte Carlo
echo ${SEP}
./monte_carlo
echo ${SEP}
echo


# Exercise 2: Page Rank
echo Exercise 2: Page Rank
echo ${SEP}
./pagerank 
echo ${SEP}
echo

# Exercise 3: Quick Sort
echo Exercise 3: Quick Sort
echo ${SEP}
./quicksort 
echo ${SEP}
echo