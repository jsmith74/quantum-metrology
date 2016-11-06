#include "MeritFunction.h"


void MeritFunction::setMeritFunction(int intParam){

    funcDimension = 2;

    return;

}



double MeritFunction::f(Eigen::VectorXd& position){

    return std::pow(position(0) * position(1) - 3,2) + 1.0;

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    std::cout << "OPTIMIZATION RESULT: " << std::endl;
    std::cout << std::pow(position(0) * position(1) - 3,2) + 1.0 << std::endl << std::endl;
    std::cout << position << std::endl << std::endl;

    return;

}



Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    return position;

}


MeritFunction::MeritFunction(){



}
