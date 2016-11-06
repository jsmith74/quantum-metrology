#include "AncillaAugment.h"

typedef Eigen::Triplet<std::complex<double> > Trip;

AncillaAugment::AncillaAugment(int N,int M,int Na,int Ma){

    photons = N;
    modes = M;

    ancillaPhotons = Na;
    ancillaModes = Ma;

    std::complex<double> Igen(0.0,1.0);
    I = Igen;

    HSdimension = g(photons,modes);
    ASdimension = g(ancillaPhotons,ancillaModes);

    initializeAugmentMatrix();

}

AncillaAugment::AncillaAugment(){

}

void AncillaAugment::setAncillaAugment(int N,int M,int Na,int Ma){

    photons = N;
    modes = M;

    ancillaPhotons = Na;
    ancillaModes = Ma;

    std::complex<double> Igen(0.0,1.0);
    I = Igen;

    HSdimension = g(photons,modes);
    ASdimension = g(ancillaPhotons,ancillaModes);

    initializeAugmentMatrix();

    return;

}


void AncillaAugment::initializeAugmentMatrix(){

    std::vector<Trip> T;

    AugmentMatrix.resize(HSdimension*ASdimension,HSdimension);

    for(int i=0;i<HSdimension*ASdimension;i++){

        T.push_back(Trip(i,i%HSdimension,1.0));

    }

    if(ASdimension == 0){

        AugmentMatrix.resize(HSdimension,HSdimension);
        for(int i=0;i<HSdimension;i++) T.push_back(Trip(i,i,1.0));

    }

    AugmentMatrix.setFromTriplets(T.begin(),T.end());

    return;

}


void AncillaAugment::setAugmentMatrix(Eigen::ArrayXd& position){

    double norm = 0.0;

    for(int i=0;i<ASdimension;i++) norm += std::pow(position(2*i),2);

    for(int i=0;i<ASdimension;i++){

        std::complex<double> ancillaValue = position(2*i) * exp(I*position(2*i+1))/sqrt(norm);
        for(int j=0;j<HSdimension;j++) AugmentMatrix.coeffRef(i*HSdimension+j,j) = ancillaValue;

    }


    return;

}


void AncillaAugment::printAugmentMatrix(){

    std::cout << "Ancilla Matrix: \n";
    std::cout << AugmentMatrix << std::endl << std::endl;
    return;
}


int AncillaAugment::g(int n,int m){
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


double AncillaAugment::doublefactorial(int x){
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
