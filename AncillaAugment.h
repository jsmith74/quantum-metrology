#ifndef ANCILLAAUGMENT_H_INCLUDED
#define ANCILLAAUGMENT_H_INCLUDED

#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Sparse>



class AncillaAugment{

    public:

        AncillaAugment(int N,int M,int Na,int Ma);
        void setAncillaAugment(int N,int M,int Na,int Ma);
        void setAugmentMatrix(Eigen::ArrayXd& position);
        void printAugmentMatrix();
        AncillaAugment();
        Eigen::SparseMatrix<std::complex<double> > AugmentMatrix;

    private:

        void initializeAugmentMatrix();
        int photons,modes,ancillaPhotons,ancillaModes,HSdimension,ASdimension;
        std::complex<double> I;
        int g(int n,int m);
        double doublefactorial(int x);

};


#endif
