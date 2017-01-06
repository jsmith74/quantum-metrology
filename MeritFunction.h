#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include "LOTransform.h"

class MeritFunction{

    public:

        MeritFunction();
        void setMeritFunction(int intParam);
        double f(Eigen::VectorXd& position);
        int funcDimension;
        void printReport(Eigen::VectorXd& position);
        Eigen::VectorXd setInitialPosition();

    private:

        Eigen::ArrayXi Row,Col;
        std::complex<double> I;
        Eigen::VectorXcd psi,psiPrime;
        Eigen::MatrixXcd OMEGAU;
        int g(int n,int m);
        double doublefactorial(int x);
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        void setRowAndCol(Eigen::ArrayXi& Row,Eigen::ArrayXi& Col,int& photons,int& modes,int& stateModes);
        int findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector,int& stateModes);
        void p_m_phiGen(Eigen::ArrayXd& p_m_phi,std::vector<Eigen::ArrayXi>& mAddress);
        void mAddressGen(std::vector<Eigen::ArrayXi>& mAddress,int& photons,int& stateModes,int& modes);
        void setmAddress(std::vector<Eigen::ArrayXi>& mAddress,Eigen::VectorXi subVector,Eigen::MatrixXi& fullVector,int& k);

};

#endif // MERITFUNCTION_H_INCLUDED
