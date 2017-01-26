#ifndef BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
#define BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED

#include "MZIMeas.h"
#include "Integration.h"
#include <iomanip>
#include <fstream>

class BranchMeasStruct{

    public:

        BranchMeasStruct();

        void setMeasChain(bool Adaptive,int numbMeas,bool import,int gridSize,double Delta);
        int setFuncDimension();
        void setPsiAndGamma(Eigen::VectorXd& position);
        void setKernalProbDistribution();

        double generalVariance();

        void printPsiAndGamma(Eigen::VectorXd& position);
        void printBranchStructure();

    private:

        Integration integrate;

        std::vector<std::vector<MZIMeas> > chainMeasurement;
        bool adaptive,import;
        int levels, numbTotalMeasBranches, numbTotalMeasOutcomes, numbGridPoints;
        double delta,dP;

        void setKernel();
        void setAdaptiveMeasurements();
        void setNonAdaptiveMeasurements();
        void setNumbTotalMeasOutcomesAdaptive();
        void setNumbTotalMeasOutcomesNonAdaptive();

};

#endif
