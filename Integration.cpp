#include "Integration.h"



Integration::Integration(){



}

void Integration::setIntegral(double Delta,double DP){

    delta = Delta;
    dP = DP;

    return;

}

double Integration::integrateArray(std::vector<double>& f){

    double output = (dP / 3.0) * f[0];

    for(int i=1;i<f.size()-2;i+=2){

        output += (4.0 * dP / 3.0) * f[i];
        output += (2.0 * dP / 3.0) * f[i+1];

    }

    output += (4.0 * dP / 3.0) * f[f.size() - 2];
    output += (dP / 3.0) * f.back();

    return output;

}
