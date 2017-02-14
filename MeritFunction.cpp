#include "MeritFunction.h"

#define PI 3.141592653589793

//#define PRINT_OUTPUT_PROB_DIST

MeritFunction::MeritFunction(){



}


Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}


void MeritFunction::setMeritFunction(int intParam,double delta,int import){

    measChain.setMeasChain(false,1,import,intParam,delta);

    measChain.setKernalProbDistribution();

    //measChain.printBranchStructure();

    funcDimension = measChain.setFuncDimension();

    deltaPrint = delta;

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    return measChain.generalVariance();

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    measChain.setPsiAndGamma(position);

    std::ofstream outfile("OptResultsDetail.dat",std::ofstream::app);

    outfile << std::endl << std::endl << std::endl;

    outfile << "OPTIMIZATION RESULT: " << std::setprecision(16) << measChain.generalVariance() << "\t" << deltaPrint << std::endl << std::endl;

    measChain.printPsiAndGamma(position,outfile);

    outfile.close();

    return;

}

void MeritFunction::printStateAmps(Eigen::VectorXd& position){

    std::ofstream outfile("BestStateAmpsDist.dat",std::ofstream::app);
    outfile << deltaPrint << "\t";
    measChain.printStateAmps(position,outfile);
    outfile << std::endl;

    outfile.close();

    outfile.open("BestGammaDist.dat",std::ofstream::app);
    outfile << deltaPrint << "\t";
    measChain.printGammaAmps(position,outfile);
    outfile << std::endl;

    outfile.close();

    measChain.setPsiAndGamma(position);

    #ifdef PRINT_OUTPUT_PROB_DIST

        measChain.printFinalProbDist();

    #endif // PRINT_OUTPUT_PROB_DIST

    return;

}

