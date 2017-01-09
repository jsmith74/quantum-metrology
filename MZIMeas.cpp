#include "MZIMeas.h"
#define PI 3.141592653589793

#define ALPHA 0.76
#define BETA 0.001

MZIMeas::MZIMeas(){


}


void MZIMeas::initializeMZIObject(int N,int M,int SM){

    std::complex<double> IGen(0.0,1.0);

    I = IGen;

    photons = N;
    modes = M;
    stateModes = SM;

    HSDimension = g(photons,modes);
    subHSDimension = g(photons,stateModes);

    psi.resize(subHSDimension);
    psiPrime.resize(HSDimension);

    setRowAndCol();
    LOOP.setLOTransform(photons,modes,Row,Col);

    mAddressGen();

    p_m_phi.resize(mAddress.size());

    branches.resize(mAddress.size());

    numbBranches = mAddress.size();

    U1 = Eigen::MatrixXcd::Identity(modes,modes);

    U23 = Eigen::MatrixXcd::Zero(modes,modes);

    UTot.resize(modes,modes);

    OMEGAU.resize(HSDimension,subHSDimension);

    initializeU23();

    return;

}

int MZIMeas::extractPhotons(){

    return photons;

}

void MZIMeas::updateP_M_PHI(){

    psiPrime = OMEGAU * psi;

    for(int i=0;i<mAddress.size();i++){

        p_m_phi(i) = norm(psiPrime(mAddress[i](0)));

        for(int j=1;j<mAddress[i].size();j++){

            p_m_phi(i) += norm(psiPrime(mAddress[i](j)));

        }

    }

    return;

}


void MZIMeas::updateOMEGAU(){

    LOOP.setUnitaryMatrixDirect(UTot);

    for(int i=0;i<OMEGAU.rows();i++){
        for(int j=0;j<OMEGAU.cols();j++){

            OMEGAU(i,j) = LOOP.omegaUij(Row(i),Col(j));

        }
    }

    return;

}


void MZIMeas::updatePhi(double& phi){

    U1(0,0) = exp(I * phi);

    UTot = U1 * U23;

    return;

}


void MZIMeas::updateGamma(double& gamma){

    U23(0,0) = cos(ALPHA) * cos(gamma);
    U23(0,1) = cos(ALPHA) * sin(gamma);
    U23(1,0) = -cos(BETA) * sin(gamma);
    U23(1,1) = cos(BETA) * cos(gamma);
    U23(2,0) = -sin(ALPHA) * cos(gamma);
    U23(2,1) = -sin(ALPHA) * sin(gamma);
    U23(3,0) = sin(BETA) * sin(gamma);
    U23(3,1) = -sin(BETA) * cos(gamma);

    return;

}

void MZIMeas::initializeU23(){

    U23(0,2) = sin(ALPHA);
    U23(1,3) = sin(BETA);
    U23(2,2) = cos(ALPHA);
    U23(3,3) = cos(BETA);

    return;

}

void MZIMeas::setPsi(Eigen::VectorXd a){

    for(int i=0;i<subHSDimension;i++){

        psi(i) = a(2*i) * exp(I * a(2*i+1));

    }

    psi.normalize();

    std::cout << "psi:\n" << std::setprecision(16)  << psi << std::endl << std::endl;
    std::cout << std::setprecision(16) << psi.norm() << std::endl;

    return;

}


int MZIMeas::findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector){

    Eigen::MatrixXi tempHold = subBasisVector.row(i);

    for(int j=0;j<fullBasisVector.rows();j++){

        for(int k=0;k<stateModes;k++){

            if(fullBasisVector(j,k) != tempHold(0,k)) break;

            if(k==stateModes-1) return j;

        }

    }

    assert(1>2 && "FAILURE");

    return -1;

}


void MZIMeas::mAddressGen(){

    int totalMeasOutComes = 0;

    for(int i=0;i<=photons;i++){

        totalMeasOutComes += g(i,stateModes);

    }

    if(stateModes == modes) totalMeasOutComes = g(photons,stateModes);

    std::cout << totalMeasOutComes << std::endl << std::endl;

    mAddress.resize(totalMeasOutComes);

    Eigen::MatrixXi fullVector = generateBasisVector(photons,modes,1);

    std::cout << "full:\n"  << fullVector << std::endl << std::endl;

    int k=0;

    for(int i=0;i<=photons;i++){

        if(stateModes == modes) i = photons;

        Eigen::MatrixXi subVector = generateBasisVector(i,stateModes,1);

        std::cout << "sub\n" << subVector << std::endl << std::endl;

        for(int j=0;j<subVector.rows();j++){

            setmAddress(subVector.row(j),fullVector,k);

            k++;

        }

    }

    return;

}

void MZIMeas::setRowAndCol(){

    Row.resize(HSDimension);

    for(int i=0;i<HSDimension;i++) Row(i) = i;

    Col.resize(subHSDimension);

    Eigen::MatrixXi fullBasisVector = generateBasisVector(photons,modes,1);
    Eigen::MatrixXi subBasisVector  = generateBasisVector(photons,stateModes,1);

    std::cout << "full\n" << fullBasisVector << std::endl << std::endl;

    std::cout << "sub start vec\n" << subBasisVector << std::endl << std::endl;

    for(int i=0;i<subHSDimension;i++){

        Col(i) = findColLoc(i,subBasisVector,fullBasisVector);

    }

    return;

}

void MZIMeas::printMathematicaMatrix(Eigen::MatrixXi& M){

    for(int i=0;i<M.rows();i++){
         std::cout << "{";
         for(int j=0;j<M.cols();j++){

            std::cout << M(i,j);
            if(j==M.cols()-1) break;
            std::cout << ",";

         }
         std::cout << "},";
    }

    return;

}

void MZIMeas::setmAddress(Eigen::VectorXi subVector,Eigen::MatrixXi& fullVector,int& k){

    for(int i=0;i<fullVector.rows();i++){

        bool match = true;

        for(int j=0;j<subVector.size();j++){

            if(subVector(j) != fullVector(i,j)){

                match = false;
                break;

            }

        }

        if(match == true){

            mAddress.at(k).conservativeResize(mAddress.at(k).size() + 1);
            mAddress.at(k)(mAddress.at(k).size() - 1) = i;

        }

    }

    return;

}


Eigen::MatrixXi MZIMeas::generateSubBasisVector(int subPhotons, int subModes){

    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    Eigen::MatrixXi nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );

    return nv;

}


Eigen::MatrixXi MZIMeas::generateBasisVector(int subPhotons,int subModes, int subMeasureModes){
    Eigen::MatrixXi output(0,subModes);
    for(int i=subPhotons;i>=0;i--){

        Eigen::MatrixXi measModes = generateSubBasisVector(i,subMeasureModes);
        if(subMeasureModes==subModes){
            return measModes;
        }

        Eigen::MatrixXi othaModes = generateSubBasisVector(subPhotons-i,subModes-subMeasureModes);

        int numbRows = measModes.rows()*othaModes.rows();
        int numbMeasRows = measModes.rows();
        int numbOthaRows = othaModes.rows();
        int outputrows = output.rows();
        output.conservativeResize(outputrows+numbRows,subModes);

        for(int k=0;k<numbMeasRows;k++){
            for(int j=0;j<numbOthaRows;j++){
                for(int l=0;l<subMeasureModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = measModes(k,l);
                }
                for(int l=subMeasureModes;l<subModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = othaModes(j,l-subMeasureModes);
                }
            }
        }
    }
    return output;
}

int MZIMeas::g(int n,int m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}


double MZIMeas::doublefactorial(int x){
    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}
