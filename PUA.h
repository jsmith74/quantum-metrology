#ifndef PUA_H_INCLUDED
#define PUA_H_INCLUDED

#include "AncillaAugment.h"
#include "LOTransform.h"

#include <fstream>

class PUA{

    public:

        PUA();
        PUA(int N,int M,int Na,int Ma,int Mm,int outcome);
        void setPUA(int N,int M,int Na,int Ma,int Mm,int outcome);
        void setQuantumOperator(Eigen::VectorXd& position);
        Eigen::MatrixXcd TotalOp;

    private:

        Eigen::MatrixXcd OMEGAU;

        int LOOPPositionSize, AAugPositionSize;

        double successProbability;

        AncillaAugment AAug;
        LOTransform LOOP;

        Eigen::ArrayXi Col,Row;
        void setCol(int N,int M,int Na,int Ma,Eigen::MatrixXi& midBasisVector,Eigen::MatrixXi& initialBasisVector);
        void setRow(Eigen::MatrixXi& midBasisVector,Eigen::MatrixXi& finalBasisVector);

        int g(int n,int m);
        double doublefactorial(int x);
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        Eigen::VectorXi measOutcome;
        void setMeasOutcome(int photons,int modes,int measModes,int outcome);
        int nonZeroCols,nonZeroRows;
        bool isEqual(Eigen::ArrayXi a1,Eigen::ArrayXi a2);

};


#endif

