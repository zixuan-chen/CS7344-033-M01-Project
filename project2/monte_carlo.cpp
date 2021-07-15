#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    int total = 10000000, count = 0;
    double x, y;
    int tn = 8;
    cout << "total trial time=" << total << endl;
    #pragma omp parallel num_threads(tn)
    {
        unsigned seed = time(NULL);

        #pragma omp for private(x, y) reduction(+:count)
        for(int i = 0; i < total; ++i){
            x = (double)rand_r(&seed) / RAND_MAX;
            y = (double)rand_r(&seed) / RAND_MAX;
            if(x*x + y*y <= 1) {
                count++;
            }
        }
    }


    double pi = 4 * (double) count / total;
    cout << "pi = " << pi << endl;
    return 0;
}
