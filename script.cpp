#include <stdlib.h>
#include <string>
#include <sstream>
#include <omp.h>
#include <unistd.h>

#define PI 3.14159265359

int main(){

    double startPoint = 0.1;
    double endPoint = PI;
    double increment = 0.05;
    int intEndPoint = (endPoint - startPoint) / increment;

    int import = 0;

#pragma omp parallel for schedule(dynamic)  default(none) \
    shared(startPoint,endPoint,increment,intEndPoint,import)
    for(int i=0;i<intEndPoint;i++){

        std::string progName = "./QuantumMetrology ";

        double delta = startPoint + i*increment;

        if( i < omp_get_num_threads() ) usleep(2000000 * omp_get_thread_num());

        std::string command,imp;
        std::stringstream ss,ssImp;
        ss << delta;
        ss >> command;
        ssImp << import;
        ssImp >> imp;

        command = progName + command + " " + imp;

        system(command.c_str());

    }

    return 0;

}
