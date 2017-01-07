#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

#include "MZIMeas.h"

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
        Eigen::VectorXcd psiPrime;
        int g(int n,int m);
        double doublefactorial(int x);

};

#endif // MERITFUNCTION_H_INCLUDED
