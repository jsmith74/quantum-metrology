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

#pragma omp parallel for schedule(dynamic)  default(none) \
    shared(startPoint,endPoint,increment,intEndPoint)
    for(int i=0;i<intEndPoint;i++){

        std::string progName = "./QuantumMetrology ";

        double delta = startPoint + i*increment;

        if( i < omp_get_num_threads() ) usleep(2000000 * omp_get_thread_num());

        std::string command;
        std::stringstream ss;
        ss << delta;
        ss >> command;
        command = progName + command;

        system(command.c_str());

    }

    return 0;

}
