#include <iostream>
#include <omp.h>
using namespace std;
#define LEN 1000000

void QuickSort(int * & array, int len){
    if(len <= 1) return;
    int pivot = array[len / 2];
    int left_ptr = 0;
    int right_ptr = len-1;
    while(left_ptr <= right_ptr){
        while(array[left_ptr] < pivot) left_ptr += 1;
        while(array[right_ptr] > pivot) right_ptr -= 1;
        if (left_ptr <= right_ptr){
            swap(array[left_ptr], array[right_ptr]);
            left_ptr += 1;
            right_ptr -= 1;
        }
    }
    int *sub_array[] = {array, &(array[left_ptr])};
    int sub_len[] = {right_ptr + 1, len - left_ptr};

#pragma omp task default(none) firstprivate(sub_array, sub_len)
    {QuickSort(sub_array[0], sub_len[0]);}
#pragma omp task default(none) firstprivate(sub_array, sub_len)
    {QuickSort(sub_array[1], sub_len[1]);}
}



int main(int argc, char const *argv[])
{
    int * array  = new int[LEN];
    clock_t startTime,endTime;
    printf("allocting a array of length %d...\n", LEN);

    if(!array){
        cout << "can't alloc enough memory\n";
    }
    for(int i = 0; i < LEN; ++i)
        array[i] = LEN - i;
    
    
    startTime = clock();
    QuickSort(array, LEN);
    endTime = clock();
    printf("Total elapsed time: %10.6f\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
    printf("checking quick sort result..\n");
    for(int i = 0; i < LEN-1; ++i){
        if (array[i] > array[i+1]){
            cout << "Quick sort incorrect" << endl;
            break;
        }
    }
    cout << "Quick sort correct\n";
    delete [] array;
    return 0;
}
