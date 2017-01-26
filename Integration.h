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

    private:

        double delta,dP;
        int numbGridPoints,levels,numbTotalMeasBranches,numbTotalMeasOutcomes;
        bool adaptive;

        std::vector<int> b,m;

        Eigen::VectorXd A,B,C;

        void setBArray(int& i,std::vector<std::vector<MZIMeas> >& chainMeasurement);
        void setMArrayAdaptive(std::vector<std::vector<MZIMeas> >& chainMeasurement);
        void printBArray();
        void printMArray();
        inline void iterateMArray(std::vector<std::vector<MZIMeas> >& chainMeasurement);

        double phi, P_phi;
        inline void updatePhiInChain(std::vector<std::vector<MZIMeas> >& chainMeasurement);
        inline void updateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double simpsonCoeff);
        inline void subUpdateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double& simpsonCoeff,int& k);

};


#endif // INTEGRATION_H_INCLUDED
