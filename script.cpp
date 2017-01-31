#include <stdlib.h>
#include <string>
#include <sstream>
#define PI 3.14159265359

int main(){

    std::string progName = "./QuantumMetrology ";

    for(double delta=0.01;delta<PI;delta+=0.02){

        std::string command;
        std::stringstream ss;
        ss << delta;
        ss >> command;
        command = progName + command;

        system(command.c_str());

    }

    return 0;

}
