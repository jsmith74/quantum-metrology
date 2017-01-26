#include "Integration.h"



Integration::Integration(){



}


void Integration::setIntegral(double Delta,double DP,int NUMBGRIDPOINTS,int LEVELS,bool ADAPTIVE,int NUMBTOTALMEASBRANCHES,int NUMBTOTALMEASOUTCOMES){

    delta = Delta;
    dP = DP;
    numbGridPoints = NUMBGRIDPOINTS;
    levels = LEVELS;
    adaptive = ADAPTIVE;
    numbTotalMeasBranches = NUMBTOTALMEASBRANCHES;
    numbTotalMeasOutcomes = NUMBTOTALMEASOUTCOMES;

    b.resize(levels);
    m.resize(levels);

    A.resize(numbTotalMeasOutcomes);
    B.resize(numbTotalMeasOutcomes);
    C.resize(numbTotalMeasBranches);

    return;

}

inline void Integration::updatePhiInChain(std::vector<std::vector<MZIMeas> >& chainMeasurement){

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement[i].size();j++){

            chainMeasurement[i][j].updatePhi(phi);
            chainMeasurement[i][j].updateOMEGAU();
            chainMeasurement[i][j].updateP_M_PHI();

        }

    }

    return;

}

inline void Integration::subUpdateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double& simpsonCoeff,int& k){

    double productHolder = P_phi;

    for(int i=0;i<levels;i++) productHolder *= chainMeasurement[i][ b[i] ].P_m_phi[ m[i] ];

    A(k) += simpsonCoeff * productHolder;

    B(k) += simpsonCoeff * phi * productHolder;

    C(k) += simpsonCoeff * phi * phi * productHolder;

    k++;

    return;

}

inline void Integration::updateABC(std::vector<std::vector<MZIMeas> >& chainMeasurement,double simpsonCoeff){

    int k=0;
    double productHolder;

    for(int i=0;i<numbTotalMeasBranches;i++){

        setBArray(i,chainMeasurement);

        if(adaptive){

            setMArrayAdaptive(chainMeasurement);

            for(int j=0;j<chainMeasurement[levels-1][i].numbBranches;j++){

                m[levels-1] = j;

                subUpdateABC(chainMeasurement,simpsonCoeff,k);

            }


        }

        else{

            for(int j=0;j<levels;j++) m[j] = 0;

            for(int j=0;j<numbTotalMeasOutcomes;j++){

                subUpdateABC(chainMeasurement,simpsonCoeff,k);

                iterateMArray(chainMeasurement);

            }

        }

    }


    return;

}

double Integration::generalVariance(std::vector<std::vector<MZIMeas> >& chainMeasurement){

    phi = -delta;
    P_phi = chainMeasurement[0][0].P_phi[0];

    updatePhiInChain(chainMeasurement);

    A = Eigen::VectorXd::Zero(numbTotalMeasOutcomes);
    B = Eigen::VectorXd::Zero(numbTotalMeasOutcomes);
    C = Eigen::VectorXd::Zero(numbTotalMeasOutcomes);           // WRITE A NEW FUNCTION THAT CIRCUMVENTS THIS INITIALIZATION (SHOULD BE REALLY EASY)

    updateABC(chainMeasurement,dP / 3.0);

    for(int i=1;i<numbGridPoints-2;i+=2){

        phi = -delta + i*dP;
        P_phi = chainMeasurement[0][0].P_phi[i];

        updatePhiInChain(chainMeasurement);

        updateABC(chainMeasurement,4.0 * dP / 3.0);


        phi = -delta + (i+1) * dP;
        P_phi = chainMeasurement[0][0].P_phi[i+1];

        updatePhiInChain(chainMeasurement);

        updateABC(chainMeasurement,2.0 * dP / 3.0);

    }

    phi = delta - dP;
    P_phi = chainMeasurement[0][0].P_phi[numbGridPoints-2];

    updatePhiInChain(chainMeasurement);

    updateABC(chainMeasurement,4.0 * dP / 3.0);

    phi = delta;
    P_phi = chainMeasurement[0][0].P_phi[numbGridPoints-1];

    updatePhiInChain(chainMeasurement);

    updateABC(chainMeasurement,dP / 3.0);

    return C.sum() - 2.0 * (B.array() / A.array()).matrix().transpose() * B + ((B.array() / A.array()) * (B.array() / A.array())).matrix().transpose() * A;

}


inline void Integration::setBArray(int& i,std::vector<std::vector<MZIMeas> >& chainMeasurement){

    b[levels-1] = i;

    for(int j=levels-1;j>0;j--){

        b[j-1] = chainMeasurement[j][ b[j] ].root;

    }

    return;

}


inline void Integration::setMArrayAdaptive(std::vector<std::vector<MZIMeas> >& chainMeasurement){

    for(int i=levels-1;i>0;i--){

        m[i-1] = chainMeasurement[i][ b[i] ].rootMeas;

    }

    return;

}

inline void Integration::iterateMArray(std::vector<std::vector<MZIMeas> >& chainMeasurement){

    m[0]++;

    for(int i=0;i<levels-1;i++){

        if(m[i] == chainMeasurement[i][0].numbBranches) {   m[i+1]++;   m[i] = 0;  }

    }

    return;

}

void Integration::printMArray(){

    std::cout << "m: ";

    for(int i=0;i<m.size();i++){

        std::cout << m[i] << " ";

    }

    std::cout << std::endl;

    return;

}

void Integration::printBArray(){

    std::cout << "b: ";

    for(int i=0;i<b.size();i++){
        std::cout << b[i] << " ";
    }

    std::cout << std::endl;

    return;

}


//double Integration::generalVariance(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m,double& phaseEstimator){
//
//    //std::ofstream test("funcIntTest.dat");
//    double phi = -delta;
//    double productHolder;
//
//    productHolder = (phi - phaseEstimator) * (phi - phaseEstimator)  * chainMeasurement[0][0].P_phi[0];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << std::setprecision(16) << 0 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    double output = (dP / 3.0) * productHolder;
//
//    for(int i=1;i<numbGridPoints-2;i+=2){
//
//        phi = -delta + i*dP;
//
//        productHolder = (phi - phaseEstimator) * (phi - phaseEstimator) * chainMeasurement[0][0].P_phi[i];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (4.0 * dP / 3.0) * productHolder;
//
//        phi = -delta + (i+1) * dP;
//
//        productHolder = (phi - phaseEstimator) * (phi - phaseEstimator) * chainMeasurement[0][0].P_phi[i+1];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i+1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (2.0 * dP / 3.0) * productHolder;
//
//    }
//
//    phi = delta - dP;
//
//    productHolder = (phi - phaseEstimator) * (phi - phaseEstimator) * chainMeasurement[0][0].P_phi[numbGridPoints-2];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-2 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (4.0 * dP / 3.0) * productHolder;
//
//
//    phi = delta;
//
//    productHolder = (phi - phaseEstimator) * (phi - phaseEstimator) * chainMeasurement[0][0].P_phi[numbGridPoints-1];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (dP / 3.0) * productHolder;
//
//    //test.close();
//
//    //std::cout << "integration result: " << std::setprecision(16) << output << std::endl;
//
//    return output;
//
//}


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


//double Integration::numer(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){
//
//    //std::ofstream test("funcIntTest.dat");
//
//    double phi = -delta;
//    double productHolder;
//
//    productHolder = phi * chainMeasurement[0][0].P_phi[0];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << std::setprecision(16) << 0 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    double output = (dP / 3.0) * productHolder;
//
//    for(int i=1;i<numbGridPoints-2;i+=2){
//
//        phi = -delta + i*dP;
//
//        productHolder = phi * chainMeasurement[0][0].P_phi[i];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (4.0 * dP / 3.0) * productHolder;
//
//        phi = -delta + (i+1) * dP;
//
//        productHolder = phi * chainMeasurement[0][0].P_phi[i+1];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i+1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (2.0 * dP / 3.0) * productHolder;
//
//    }
//
//    phi = delta - dP;
//
//    productHolder = phi * chainMeasurement[0][0].P_phi[numbGridPoints-2];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-2 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (4.0 * dP / 3.0) * productHolder;
//
//
//    phi = delta;
//
//    productHolder = phi * chainMeasurement[0][0].P_phi[numbGridPoints-1];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (dP / 3.0) * productHolder;
//
//    //test.close();
//
//    //std::cout << "integration result: " << std::setprecision(16) << output << std::endl;
//
//    return output;
//
//}
//
//
//double Integration::denom(std::vector<std::vector<MZIMeas> >& chainMeasurement,std::vector<int>& b,std::vector<int>& m){
//
//    //std::ofstream test("funcIntTest.dat");
//
//    double phi = -delta;
//    double productHolder;
//
//    productHolder = chainMeasurement[0][0].P_phi[0];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << std::setprecision(16) << 0 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    double output = (dP / 3.0) * productHolder;
//
//    for(int i=1;i<numbGridPoints-2;i+=2){
//
//        phi = -delta + i*dP;
//
//        productHolder = chainMeasurement[0][0].P_phi[i];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (4.0 * dP / 3.0) * productHolder;
//
//        phi = -delta + (i+1) * dP;
//
//        productHolder = chainMeasurement[0][0].P_phi[i+1];
//
//        for(int j=0;j<levels;j++){
//
//            chainMeasurement[j][b[j]].updatePhi(phi);
//            chainMeasurement[j][b[j]].updateOMEGAU();
//            chainMeasurement[j][b[j]].updateP_M_PHI(m[j]);
//
//            productHolder *= chainMeasurement[j][b[j]].P_m_phi[m[j]];
//
//        }
//
//        //test << i+1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//        output += (2.0 * dP / 3.0) * productHolder;
//
//    }
//
//    phi = delta - dP;
//
//    productHolder = chainMeasurement[0][0].P_phi[numbGridPoints-2];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-2 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (4.0 * dP / 3.0) * productHolder;
//
//
//    phi = delta;
//
//    productHolder = chainMeasurement[0][0].P_phi[numbGridPoints-1];
//
//    for(int i=0;i<levels;i++){
//
//        chainMeasurement[i][b[i]].updatePhi(phi);
//        chainMeasurement[i][b[i]].updateOMEGAU();
//        chainMeasurement[i][b[i]].updateP_M_PHI(m[i]);
//
//        productHolder *= chainMeasurement[i][b[i]].P_m_phi[m[i]];
//
//    }
//
//    //test << numbGridPoints-1 << "\t" << phi << "\t" << productHolder << std::endl;
//
//    output += (dP / 3.0) * productHolder;
//
//    //test.close();
//
//    //std::cout << "integration result: " << std::setprecision(16) << output << std::endl;
//
//    //if(m[0] == 4 && m[1] == 4 && m[2] == 4) assert(1>2 && "UP TO HERE");
//
//    // UP TO HERE - CHECK THAT THE INTEGRATION IS WORKING PROPERLY - GO THROUGH THIS FUNCTION AN EXTRA TIME JUST TO MAKE SURE EVERYTHING IS GOOD
//    // SOMETHING ISN'T RIGHT - CHECK IT
//
//
//    return output;
//
//}



