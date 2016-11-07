#include "MeritFunction.h"

#define PI 3.141592653589793


void MeritFunction::setMeritFunction(int intParam){

    photons = 2;
    realModes = 2;
    lossModes = 2;
    interferometerControlParams = 1;

    modes = realModes + lossModes;

    std::complex<double> IGen(0.0,1.0);

    I = IGen;

    U1 = Eigen::MatrixXcd::Identity(modes,modes);

    U1(0,0) = std::exp(I * PI/4.0);

    double alpha,beta;
    alpha = 0.0;
    beta = 0.0;
    U2.resize(modes,modes);

    U2 << cos(alpha),0.0,sin(alpha),0.0,
            0.0,cos(beta),0.0,sin(beta),
            -sin(alpha), 0.0,cos(alpha),0.0,
            0.0,-sin(beta),0.0,cos(beta);

    U3 = Eigen::MatrixXcd::Identity(modes,modes);

    HSDimension = g(photons,modes);
    subHSDimension = g(photons,realModes);

    funcDimension = interferometerControlParams + 2 * subHSDimension;

    initialState = Eigen::VectorXcd::Zero(subHSDimension);

    setRowAndCol(Row,Col);

    LOOP.setLOTransform(photons,modes,Row,Col);

    OMEGAU.resize(HSDimension,subHSDimension);

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    for(int i=0;i<subHSDimension;i++) initialState(i) = position(2*i) * exp(I * position(2*i+1));

    initialState.normalize();

    double gamma = position(2*subHSDimension);

    U3(0,0) = cos(gamma);       U3(0,1) = sin(gamma);
    U3(1,0) = -sin(gamma);      U3(1,1) = cos(gamma);

    UTot = U1 * U2 * U3;

    LOOP.setUnitaryMatrixDirect(UTot);

    for(int i=0;i<HSDimension;i++){
        for(int j=0;j<subHSDimension;j++){

            OMEGAU(i,j) = LOOP.omegaUij(Row(i),Col(j));

        }
    }

    finalState = OMEGAU * initialState;

    std::cout << initialState << std::endl << std::endl;

    std::cout << finalState << std::endl << std::endl;

    /** ==== TO DO -- CHECK THAT THIS TRANSFORMATION MATCHES YANG's ==== */

    assert(1>2 && "break it here");

    return std::pow(position(0) * position(1) - 3,2) + 1.0;

}


void MeritFunction::setRowAndCol(Eigen::ArrayXi& Row,Eigen::ArrayXi& Col){

    Row.resize(HSDimension);

    for(int i=0;i<HSDimension;i++) Row(i) = i;

    Col.resize(subHSDimension);

    Eigen::MatrixXi fullBasisVector = generateBasisVector(photons,modes,1);
    Eigen::MatrixXi subBasisVector  = generateBasisVector(photons,realModes,1);

    for(int i=0;i<subHSDimension;i++){

        Col(i) = findColLoc(i,subBasisVector,fullBasisVector);

    }

    return;

}

int MeritFunction::findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector){

    Eigen::MatrixXi tempHold = subBasisVector.row(i);

    for(int j=0;j<fullBasisVector.rows();j++){

        for(int k=0;k<realModes;k++){

            if(fullBasisVector(j,k) != tempHold(0,k)) break;

            if(k==realModes-1) return j;

        }

    }

    assert(1>2 && "FAILURE");

    return -1;

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    std::cout << "OPTIMIZATION RESULT: " << std::endl;
    std::cout << std::pow(position(0) * position(1) - 3,2) + 1.0 << std::endl << std::endl;
    std::cout << position << std::endl << std::endl;

    return;

}



Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}


MeritFunction::MeritFunction(){



}

Eigen::MatrixXi MeritFunction::generateBasisVector(int subPhotons,int subModes, int subMeasureModes){
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

int MeritFunction::g(int n,int m){
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


double MeritFunction::doublefactorial(int x){
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



Eigen::MatrixXi MeritFunction::generateSubBasisVector(int subPhotons, int subModes){
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
