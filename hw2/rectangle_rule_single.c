#include <stdio.h>
/* This program computes pi using the rectangle rule. */ 
#define INTERVALS 1000000 
int main (int argc, char *argv[]) 
{ 
    double area; /* Area under curve */ 
    double ysum; /* Sum of rectangle heights */ 
    double xi; /* Midpoint of interval */ 
    int i; 
    ysum = 0.0; 
    for (i = 0; i < INTERVALS; i++) { 
        xi = (1.0/INTERVALS)*(i+0.5); 
        ysum += 4.0/(1.0+xi*xi); 
    } 
    area = ysum * (1.0 / INTERVALS); 
    printf ("Area is %13.11f\n", area); 
    return 0; 
} 
