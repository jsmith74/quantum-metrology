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
        void printBranchStructure();
        int setFuncDimension();
        void setPsiAndGamma(Eigen::VectorXd& position);
        void setPhaseEstimators();
        void setKernalProbDistribution();

    private:

        std::vector<int> m;
        std::vector<int> b;
        std::vector<double> phaseEstimators;
        std::vector<std::vector<MZIMeas> > chainMeasurement;
        bool adaptive,import;
        int levels, numbTotalMeasBranches, numbTotalMeasOutcomes, numbGridPoints;
        double delta,dP;

        void setKernel();
        void setAdaptiveMeasurements();
        void setNonAdaptiveMeasurements();
        void setNumbTotalMeasOutcomesAdaptive();
        void setNumbTotalMeasOutcomesNonAdaptive();
        inline void setBArray(int& i);
        void printBArray();
        inline void setMArrayAdaptive(int& j);
        void printMArray();

};

#endif // BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
