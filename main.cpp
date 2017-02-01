#include "BFGS_Optimization.h"

#include <time.h>
#include <fstream>
#include <iostream>

int main( int argc, char *argv[] ){

    if(argc != 2){

        std::cout << "./QuantumMetrology [uncertainty inverval size]" << std::endl;
        return 1;

    }

    double gradientCheck = 1e-8;

    double maxStepSize = 2.0;

    int integrationGridSize = 6000;

    int optimizationAttempts = 10;

    double bestResult = 1e20;

    double delta = std::atof(argv[1]);

    clock_t t1,t2;

    t1 = clock();

    for(int i=0;i<optimizationAttempts;i++){

        BFGS_Optimization optimizer(gradientCheck,maxStepSize,integrationGridSize,delta);

        double result = optimizer.minimize();

        std::ofstream outfile("OptResults.dat",std::ofstream::app);

        outfile << i << "\t" << delta << "\t" << std::setprecision(16) << result << std::endl;

        outfile.close();

        if(result < bestResult) {

            bestResult = result;
            optimizer.printStateAmps();  // TO DO: FIX THIS SO IT ONLY PRINTS THE TRULY BEST ONE

        }

    }

    t2 = clock();

    float diff = (float)t2 - (float)t1;

    //std::cout << "Runtime: " << diff/CLOCKS_PER_SEC << std::endl << std::endl;

    std::ofstream outfile("BestResults.dat",std::ofstream::app);
    outfile << delta << "\t" << std::setprecision(16) << bestResult << std::endl;
    outfile.close();

    return 0;

}
