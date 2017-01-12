#ifndef BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
#define BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED

#include "MZIMeas.h"
#include <fstream>

class BranchMeasStruct{

    public:

        BranchMeasStruct();
        void setMeasChain(bool Adaptive,int numbMeas,bool import);
        void printBranchStructure();
        int setFuncDimension();
        void setPsiAndGamma(Eigen::VectorXd& position);
        void setPhaseEstimators();

    private:

        std::vector<int> m;
        std::vector<double> phaseEstimators;
        std::vector<std::vector<MZIMeas> > chainMeasurement;
        bool adaptive,import;
        int levels, numbTotalMeasBranches, numbTotalMeasOutcomes;

        void setKernel();
        void setAdaptiveMeasurements();
        void setNonAdaptiveMeasurements();
        void setNumbTotalMeasOutcomesAdaptive();
        void setNumbTotalMeasOutcomesNonAdaptive();

};

#endif // BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
