#ifndef MZIMEAS_H_INCLUDED
#define MZIMEAS_H_INCLUDED

#include <Eigen/Dense>
#include <vector>
#include <iomanip>

#include "LOTransform.h"

class MZIMeas{

    public:

        MZIMeas();

        void initializeMZIObject(int N,int M,int SM);
        void setPsi(Eigen::VectorXd a);
        void updateGamma(double& gamma);
        void updatePhi(double& phi);
        void updateOMEGAU();
        void updateP_M_PHI();

        Eigen::ArrayXd p_m_phi;

    private:

        std::complex<double> I;
        int photons,modes,stateModes,HSDimension,subHSDimension;
        double phi,gamma;
        Eigen::VectorXcd psi,psiPrime;
        Eigen::MatrixXcd OMEGAU,UTot,U1,U23;
        Eigen::ArrayXi Row,Col;
        std::vector<Eigen::ArrayXi> mAddress;
        LOTransform LOOP;  //REWRITE THIS SO THAT INTERFEROMETERS WITH IDENTICAL STRUCTURE SHARE LOOTRANSFORM OBJECTS (IF THIS TURNS OUT AS BOTTLENECK)

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

#endif // MZIMEAS_H_INCLUDED
