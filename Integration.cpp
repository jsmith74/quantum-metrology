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

double Integration::generalVariance(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<double>& phaseEstimators){

    return 2.0;

}


double Integration::numer(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){

    //std::ofstream test("funcIntTest.dat");

    double phi = -delta;
    double productHolder;

    productHolder = phi * chainMeasurement[0][0].P_phi[0];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << std::setprecision(16) << 0 << "\t" << phi << "\t" << productHolder << std::endl;

    double output = (dP / 3.0) * productHolder;

    for(int i=1;i<numbGridPoints-2;i+=2){

        phi = -delta + i*dP;

        productHolder = phi * chainMeasurement[0][0].P_phi[i];

        for(int j=0;j<levels;j++){

            chainMeasurement[j][b[j]].updatePhi(phi);
            chainMeasurement[j][b[j]].updateOMEGAU();
            chainMeasurement[j][b[j]].updateP_M_PHI();

            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];

        }

        //test << i << "\t" << phi << "\t" << productHolder << std::endl;

        output += (4.0 * dP / 3.0) * productHolder;

        phi = -delta + (i+1) * dP;

        productHolder = phi * chainMeasurement[0][0].P_phi[i+1];

        for(int j=0;j<levels;j++){

            chainMeasurement[j][b[j]].updatePhi(phi);
            chainMeasurement[j][b[j]].updateOMEGAU();
            chainMeasurement[j][b[j]].updateP_M_PHI();

            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];

        }

        //test << i+1 << "\t" << phi << "\t" << productHolder << std::endl;

        output += (2.0 * dP / 3.0) * productHolder;

    }

    phi = delta - dP;

    productHolder = phi * chainMeasurement[0][0].P_phi[numbGridPoints-2];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << numbGridPoints-2 << "\t" << phi << "\t" << productHolder << std::endl;

    output += (4.0 * dP / 3.0) * productHolder;


    phi = delta;

    productHolder = phi * chainMeasurement[0][0].P_phi[numbGridPoints-1];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << numbGridPoints-1 << "\t" << phi << "\t" << productHolder << std::endl;

    output += (dP / 3.0) * productHolder;

    //test.close();

    //std::cout << "integration result: " << std::setprecision(16) << output << std::endl;

    return output;

}


double Integration::denom(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){

    //std::ofstream test("funcIntTest.dat");

    double phi = -delta;
    double productHolder;

    productHolder = chainMeasurement[0][0].P_phi[0];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << std::setprecision(16) << 0 << "\t" << phi << "\t" << productHolder << std::endl;

    double output = (dP / 3.0) * productHolder;

    for(int i=1;i<numbGridPoints-2;i+=2){

        phi = -delta + i*dP;

        productHolder = chainMeasurement[0][0].P_phi[i];

        for(int j=0;j<levels;j++){

            chainMeasurement[j][b[j]].updatePhi(phi);
            chainMeasurement[j][b[j]].updateOMEGAU();
            chainMeasurement[j][b[j]].updateP_M_PHI();

            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];

        }

        //test << i << "\t" << phi << "\t" << productHolder << std::endl;

        output += (4.0 * dP / 3.0) * productHolder;

        phi = -delta + (i+1) * dP;

        productHolder = chainMeasurement[0][0].P_phi[i+1];

        for(int j=0;j<levels;j++){

            chainMeasurement[j][b[j]].updatePhi(phi);
            chainMeasurement[j][b[j]].updateOMEGAU();
            chainMeasurement[j][b[j]].updateP_M_PHI();

            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];

        }

        //test << i+1 << "\t" << phi << "\t" << productHolder << std::endl;

        output += (2.0 * dP / 3.0) * productHolder;

    }

    phi = delta - dP;

    productHolder = chainMeasurement[0][0].P_phi[numbGridPoints-2];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << numbGridPoints-2 << "\t" << phi << "\t" << productHolder << std::endl;

    output += (4.0 * dP / 3.0) * productHolder;


    phi = delta;

    productHolder = chainMeasurement[0][0].P_phi[numbGridPoints-1];

    for(int i=0;i<levels;i++){

        chainMeasurement[i][b[i]].updatePhi(phi);
        chainMeasurement[i][b[i]].updateOMEGAU();
        chainMeasurement[i][b[i]].updateP_M_PHI();

        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];

    }

    //test << numbGridPoints-1 << "\t" << phi << "\t" << productHolder << std::endl;

    output += (dP / 3.0) * productHolder;

    //test.close();

    //std::cout << "integration result: " << std::setprecision(16) << output << std::endl;

    //if(m[0] == 4 && m[1] == 4 && m[2] == 4) assert(1>2 && "UP TO HERE");

    // UP TO HERE - CHECK THAT THE INTEGRATION IS WORKING PROPERLY - GO THROUGH THIS FUNCTION AN EXTRA TIME JUST TO MAKE SURE EVERYTHING IS GOOD
    // SOMETHING ISN'T RIGHT - CHECK IT


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
