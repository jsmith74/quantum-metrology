#ifndef BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
#define BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED

#include "MZIMeas.h"
#include "Integration.h"
#include <iomanip>
#include <fstream>

class BranchMeasStruct{

    public:

        BranchMeasStruct();

        void setMeasChain(bool Adaptive,int numbMeas,int Import,int gridSize,double Delta);
        int setFuncDimension();
        void setPsiAndGamma(Eigen::VectorXd& position);
        void setKernalProbDistribution();

        double generalVariance();

        void printPsiAndGamma(Eigen::VectorXd& position,std::ofstream& outfile);
        void printBranchStructure();

        void printStateAmps(Eigen::VectorXd& position,std::ofstream& outfile);
        void printGammaAmps(Eigen::VectorXd& position,std::ofstream& outfile);
        void printRelativePhase(Eigen::VectorXd& position,std::ofstream& outfile);

        void printFinalProbDist();

    private:

        Integration integrate;

        std::vector<std::vector<MZIMeas> > chainMeasurement;
        bool adaptive;

        int levels, numbTotalMeasBranches, numbTotalMeasOutcomes, numbGridPoints, import;
        double delta,dP;

        void setKernel();
        void setAdaptiveMeasurements();
        void setNonAdaptiveMeasurements();
        void setNumbTotalMeasOutcomesAdaptive();
        void setNumbTotalMeasOutcomesNonAdaptive();

        void setNumbPrevBranches(unsigned long& numbPrevBranches,std::ifstream& infile);

};

#endif
