#include "MeritFunction.h"

#define PI 3.141592653589793
#define ALPHA 0.0;
#define BETA 0.0

Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd position = Eigen::VectorXd::Random(funcDimension);

    position *= PI;

    return position;

}

void MeritFunction::setMeritFunction(int intParam){

    std::complex<double> IGen(0.0,1.0);

    I = IGen;

    /** =========================================== */

    funcDimension = 3;

    int photons = 2;

    int stateModes = 2;
    int modes = 4;

    setRowAndCol(Row,Col,photons,modes,stateModes);

    std::cout << "col\n" << Col << std::endl << std::endl;

    LOTransform LOOPtest;
    LOOPtest.setLOTransform(photons,modes,Row,Col);

    std::vector<Eigen::ArrayXi> mAddressTest;
    mAddressGen(mAddressTest,photons,stateModes,modes);

    Eigen::MatrixXi fullVector = generateBasisVector(photons,modes,1);

    for(int i=0;i<mAddressTest.size();i++){

        std::cout << i << ":\n";
        for(int j=0;j<mAddressTest.at(i).size();j++) std::cout << mAddressTest.at(i)(j) << " ";
        std::cout << std::endl;
        for(int j=0;j<mAddressTest.at(i).size();j++) std::cout << fullVector.row(mAddressTest.at(i)(j))  << std::endl;
        std::cout << std::endl;
    }

    // UP TO HERE:
    // CHECK THAT THIS CODE IS PRODUCING THE CORRECT MAPS FROM STARTING STATES ONTO FULL HILBERT SPACE
    // AND THE CORRECT MEASUREMENT LOCATION MAPS FOR RANDOM PHYSICAL SYSTEMS (MODES, PHOTONS, STATEMODES)

    /** =========================================== */

    assert(1>2 && "End here");

    return;

}

void MeritFunction::mAddressGen(std::vector<Eigen::ArrayXi>& mAddress,int& photons,int& stateModes,int& modes){

    int totalMeasOutComes = 0;

    for(int i=0;i<=photons;i++){

        totalMeasOutComes += g(i,stateModes);

    }

    std::cout << totalMeasOutComes << std::endl << std::endl;

    mAddress.resize(totalMeasOutComes);

    Eigen::MatrixXi fullVector = generateBasisVector(photons,modes,1);

    std::cout << "full:\n"  << fullVector << std::endl << std::endl;

    int k=0;

    for(int i=0;i<=photons;i++){

        Eigen::MatrixXi subVector = generateBasisVector(i,stateModes,1);

        std::cout << "sub\n" << subVector << std::endl << std::endl;

        for(int j=0;j<subVector.rows();j++){

            setmAddress(mAddress,subVector.row(j),fullVector,k);

            k++;

        }

    }

    return;

}

void MeritFunction::setmAddress(std::vector<Eigen::ArrayXi>& mAddress,Eigen::VectorXi subVector,Eigen::MatrixXi& fullVector,int& k){

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

void MeritFunction::p_m_phiGen(Eigen::ArrayXd& p_m_phi,double& phi,double& gamma,Eigen::VectorXcd& psi,LOTransform& LOOP){

    // TURN THIS INTO A FUNCTION AS A PART OF A SEPARATE OBJECT WITH mAddress and all this stuff as data. I THINK
    // EACH p(m|phi) WILL BE A SEPARATE OBJECT OF THE SAME TYPE REPRESENTING A DIFFERENT RUN THROUGH INTERFEROMETER

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    return 1.0;

}


void MeritFunction::setRowAndCol(Eigen::ArrayXi& Row,Eigen::ArrayXi& Col,int& photons,int& modes,int& stateModes){

    int HSDimension = g(photons,modes);

    int subHSDimension = g(photons,stateModes);

    Row.resize(HSDimension);

    for(int i=0;i<HSDimension;i++) Row(i) = i;

    Col.resize(subHSDimension);

    Eigen::MatrixXi fullBasisVector = generateBasisVector(photons,modes,1);
    Eigen::MatrixXi subBasisVector  = generateBasisVector(photons,stateModes,1);

    std::cout << "full\n" << fullBasisVector << std::endl << std::endl;
    std::cout << "sub start vec\n" << subBasisVector << std::endl << std::endl;

    for(int i=0;i<subHSDimension;i++){

        Col(i) = findColLoc(i,subBasisVector,fullBasisVector,stateModes);

    }

    return;

}

int MeritFunction::findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector,int& stateModes){

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


void MeritFunction::printReport(Eigen::VectorXd& position){

    std::cout << "OPTIMIZATION RESULT: " << std::endl;
    std::cout << std::pow(position(0) * position(1) - 3,2) + 1.0 << std::endl << std::endl;
    std::cout << position << std::endl << std::endl;

    return;

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
