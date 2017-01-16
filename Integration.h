#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include <vector>
#include <iostream>

class Integration{

    public:

        Integration();
        void setIntegral(double Delta,double DP);
        double integrateArray(std::vector<double>& f);

    private:

        double delta,dP;

};


#endif // INTEGRATION_H_INCLUDED
