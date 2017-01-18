#include "Integration.h"



Integration::Integration(){



}


void Integration::setIntegral(double Delta,double DP,int NUMBGRIDPOINTS,int LEVELS){

    delta = Delta;
    dP = DP;
    numbGridPoints = NUMBGRIDPOINTS;
    levels = LEVELS;

    return;

}


double Integration::numer(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){

    double output;

    return output;

}


double Integration::denom(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){

    double phi = -delta;
    double productHolder = 1.0;

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    double output = (dP / 3.0) * ;

    return output;

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
