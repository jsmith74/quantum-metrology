#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

#include "MZIMeas.h"
#include "BranchMeasStruct.h"

class MeritFunction{

    public:

        MeritFunction();

        Eigen::VectorXd setInitialPosition();
        void setMeritFunction(int intParam,double delta,int import);
        double f(Eigen::VectorXd& position);
        void printReport(Eigen::VectorXd& position);

        int funcDimension;

        void printStateAmps(Eigen::VectorXd& position);

    private:

        Eigen::ArrayXi Row,Col;
        Eigen::VectorXcd psiPrime;
        BranchMeasStruct measChain;
        double deltaPrint;

};

#endif
