#ifndef MZIMEAS_H_INCLUDED
#define MZIMEAS_H_INCLUDED

#include <Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>

#include "LOTransform.h"

class MZIMeas{

    public:

        MZIMeas();

        void setPsi(Eigen::VectorXd& position,int& k);
        void updateGamma(double& gamma);
        void updatePhi(double& phi);
        void updateOMEGAU();
        void updateP_M_PHI();
        void initializeMZIObject(int N,int M,int SM);

        double delta,dP;
        std::vector<double> P_m_phi;
        std::vector<double> P_phi;
        int photons,modes,stateModes,HSDimension,subHSDimension,level,root,numbBranches,rootMeas;

        void printPDist();
        void printMAddress();
        void printPsi(Eigen::VectorXd& position,int& k,std::ofstream& outfile);
        int extractPhotons();

    private:

        std::complex<double> I;
        Eigen::VectorXcd psi,psiPrime;
        Eigen::MatrixXcd OMEGAU,UTot,U1,U23;
        Eigen::ArrayXi Row,Col;
        std::vector<Eigen::ArrayXi> mAddress;
        LOTransform LOOP;

        int g(int n,int m);
        double doublefactorial(int x);
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        void setRowAndCol();
        void mAddressGen();
        void setmAddress(Eigen::VectorXi subVector,Eigen::MatrixXi& fullVector,int& k);
        int findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector);
        void initializeU23();

};

#endif
