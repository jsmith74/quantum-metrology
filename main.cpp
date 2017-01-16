#include "BFGS_Optimization.h"

#include <time.h>
#include <fstream>
#include <iostream>


int main( int argc, char *argv[] ){

    double gradientCheck = 1e-6;

    double maxStepSize = 2.0;

    int intParam = 6000;

    clock_t t1,t2;

    BFGS_Optimization optimizer(gradientCheck,maxStepSize,intParam);

    t1 = clock();
    optimizer.minimize();
    t2 = clock();

    float diff = (float)t2 - (float)t1;

    std::cout << "Runtime: " << diff/CLOCKS_PER_SEC << std::endl << std::endl;

    return 0;

}
