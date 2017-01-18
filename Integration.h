#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include <vector>
#include <iostream>
#include "MZIMeas.h"

class Integration{

    public:

        Integration();
        void setIntegral(double Delta,double DP,int S,int LEVELS);
        double integrateArray(std::vector<double>& f);
        double numer(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m);
        double denom(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m);

    private:

        double delta,dP;
        int numbGridPoints,levels;

};


#endif // INTEGRATION_H_INCLUDED
