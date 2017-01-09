#include "MeritFunction.h"

#define PI 3.141592653589793

Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}

void MeritFunction::setMeritFunction(int intParam){

    funcDimension = 3;

    measChain.setMeasChain(true,3,false);

    measChain.printBranchStructure();

    // UP TO HERE PLAY WITH BRANCH STRUCTURE, MAKE SURE ITS GENERATING THE RIGHT STRUCTURE AND MAKE SURE MEASUREMENT OUTCOMES ARE WORKING
    // DONT FORGET TO CHECK FOR NON-ADAPTIVE REPEATED MEASUREMENTS
    // IF ALL IS GOOD, SET ALL N=2 FOR ALL BRANCHES AND MAYBE WRITE AN I/O INTERFACE FOR INPUTTING N
    // THEN, TRY AND MEET WITH KAPLAN TO FIND OUT WHAT HAPPENS TO P(PHI) AFTER SUCCESSIVE MEASUREMENTS

//    MZIMeas MZITest;
//
//    MZITest.initializeMZIObject(2,4,2);
//
//    MZITest.setPsi(Eigen::VectorXd::Random(6));
//
//    double phiTest = PI / 6.0;
//
//    double gammaTest = 0.54;
//
//    MZITest.updateGamma(gammaTest);
//
//    MZITest.updatePhi(phiTest);
//
//    MZITest.updateOMEGAU();
//
//    MZITest.updateP_M_PHI();

    assert(1>2 && "End here");

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    return 1.0;

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    std::cout << "OPTIMIZATION RESULT: " << std::endl;
    std::cout << std::pow(position(0) * position(1) - 3,2) + 1.0 << std::endl << std::endl;
    std::cout << position << std::endl << std::endl;

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


