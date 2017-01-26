#include "MeritFunction.h"

#define PI 3.141592653589793


MeritFunction::MeritFunction(){



}


Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}


void MeritFunction::setMeritFunction(int intParam){

    measChain.setMeasChain(true,3,false,intParam,.1);

    measChain.setKernalProbDistribution();

    measChain.printBranchStructure();

    funcDimension = measChain.setFuncDimension();

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    return measChain.generalVariance();

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    std::cout << "OPTIMIZATION RESULT: " << measChain.generalVariance() << std::endl << std::endl;

    measChain.printPsiAndGamma(position);

    return;

}
