#include "BFGS_Optimization.h"

#include <time.h>
#include <fstream>
#include <iostream>
#include <omp.h>

int main( int argc, char *argv[] ){

    double gradientCheck = 1e-8;

    double maxStepSize = 2.0;

    int integrationGridSize = 6000;

    int optimizationAttempts = 100;

    clock_t t1,t2;

    t1 = clock();

#pragma omp parallel for schedule(dynamic)  default(none) \
    shared(gradientCheck,maxStepSize,integrationGridSize,optimizationAttempts)
    for(int i=0;i<optimizationAttempts;i++){

        BFGS_Optimization optimizer(gradientCheck,maxStepSize,integrationGridSize);

        double result = optimizer.minimize();

        std::ofstream outfile("OptResults.dat",std::ofstream::app);

        outfile << i << "\t" << omp_get_thread_num() << "\t" << result << std::endl;

        outfile.close();

    }

    t2 = clock();

    float diff = (float)t2 - (float)t1;

    std::cout << "Runtime: " << diff/CLOCKS_PER_SEC << std::endl << std::endl;

    return 0;

}
