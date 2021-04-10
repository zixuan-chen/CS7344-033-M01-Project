#include <iostream>
#include <stdio.h>
using namespace std;

int main(){
    unsigned char a[5] = {13, 22, 43, 64, 99};
    unsigned char add = 0, multiply = 1, maximum = 0, minimum = 255, bitwise_or = 0, bitwise_and = 1, logicalOR = 0, logicalAND = 1;
    for (int i = 0; i < 5; i++)
    {
        add += a[i];
        multiply *= a[i];
        maximum = max(maximum, a[i]);
        minimum = min(minimum, a[i]);
        bitwise_or = bitwise_or | a[i];
        bitwise_and = bitwise_and & a[i];
        logicalOR = logicalOR || a[i];
        logicalAND = logicalAND && a[i];

    }

    printf("add = %u\n multiply = %u\n maximum = %u\n minimum = %u\n bitwise_or = %u\n bitwise_and = %u\n logicalOR = %u\n logicalAND = %u\n", 
            add, multiply, maximum, minimum, bitwise_or, bitwise_and, logicalOR, logicalAND);
    return 0;
}