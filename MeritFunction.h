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
        Eigen::MatrixXcd UTot;
        std::complex<double> I;
        Eigen::MatrixXcd OMEGAU;
        int HSDimension,subHSDimension,photons,interferometerControlParams,modes,stateModes;
        Eigen::VectorXcd initialState,finalState;
        int g(int n,int m);
        double doublefactorial(int x);
        LOTransform LOOP;
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        void setRowAndCol(Eigen::ArrayXi& Row,Eigen::ArrayXi& Col);
        int findColLoc(int i,Eigen::MatrixXi& subBasisVector,Eigen::MatrixXi& fullBasisVector);
};

#endif // MERITFUNCTION_H_INCLUDED
