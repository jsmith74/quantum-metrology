#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include <vector>
#include <iostream>
#include <iomanip>
#include "MZIMeas.h"

class Integration{

    public:

        Integration();

        void setIntegral(double Delta,double DP,int S,int LEVELS,bool ADAPTIVE,int NUMBTOTALMEASBRANCHES,int NUMBTOTALMEASOUTCOMES);

        double integrateArray(std::vector<double>& f);
        double generalVariance(std::vector<std::vector<MZIMeas> >& chainMeasurement);

        void printFinalProbDist(std::vector<std::vector<MZIMeas> >& chainMeasurement);

    private:

        bool adaptive;
        double delta,dP,phi, P_phi;
        int numbGridPoints,levels,numbTotalMeasBranches,numbTotalMeasOutcomes;

        std::string filename;

        std::vector<int> b,m;

        Eigen::VectorXd A,B,C;

        inline void iterateMArray(std::vector<std::vector<MZIMeas> >& chainMeasurement);
        inline void setBArray(int& i,std::vector<std::vector<MZIMeas> >& chainMeasurement);
        inline void setMArrayAdaptive(std::vector<std::vector<MZIMeas> >& chainMeasurement);
        inline void updatePhiInChain(std::vector<std::vector<MZIMeas> >& chainMeasurement);
        inline void updateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double simpsonCoeff);
        inline void subUpdateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double& simpsonCoeff,int& k);
        inline void initializeABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double simpsonCoeff);
        inline void subInitializeABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double& simpsonCoeff,int& k);

        void printBArray();
        void printMArray();

        inline void subPrintFinalProbDist(std::vector<std::vector<MZIMeas> >& chainMeasurement);

};


#endif
