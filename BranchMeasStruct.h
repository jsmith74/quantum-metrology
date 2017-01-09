#ifndef BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
#define BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED

#include "MZIMeas.h"
#include <fstream>

class BranchMeasStruct{

    public:

        BranchMeasStruct();
        void setMeasChain(bool Adaptive,int numbMeas,bool import);
        void printBranchStructure();

    private:

        std::vector<std::vector<MZIMeas> > chainMeasurement;
        bool adaptive,import;
        int levels;

        void setKernel();
        void setAdaptiveMeasurements();
        void setNonAdaptiveMeasurements();

};

#endif // BRANCHMEASUREMENTSTRUCTURE_H_INCLUDED
