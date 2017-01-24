#include "MeritFunction.h"

#define PI 3.141592653589793

Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}

void MeritFunction::setMeritFunction(int intParam){

    measChain.setMeasChain(true,1,false,intParam,.1);

    measChain.setKernalProbDistribution();

    measChain.printBranchStructure();

    funcDimension = measChain.setFuncDimension();

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    measChain.setPhaseEstimators();

    return measChain.generalVariance();

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    measChain.setPhaseEstimators();

    std::cout << "OPTIMIZATION RESULT: " << measChain.generalVariance() << std::endl << std::endl;

    measChain.printPsiAndGamma(position);

    return;

}



MeritFunction::MeritFunction(){



}



int MeritFunction::g(int n,int m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}


double MeritFunction::doublefactorial(int x){
    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}


