#ifndef LOTRANSFORM_H_INCLUDED
#define LOTRANSFORM_H_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include <iostream>


class LOTransform{

    public:

        std::complex<double> omegaUij(int& i,int& j);
        LOTransform(int N,int M,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col);
        void setLOTransform(int N,int M,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col);
        void setUnitaryMatrix(Eigen::ArrayXd& position);
        LOTransform();

    private:

        int photons, modes;
        std::complex<double> I;
        Eigen::Array2i LDimension;
        Eigen::MatrixXcd U;
        std::vector<Eigen::ArrayXd> CijMatrix;
        std::vector<Eigen::ArrayXi> kMatrix;
        inline int genKIndex(int& i,int& j,int& dim);
        inline std::complex<double> UProduct(int& s,int& p,int kIndex);
        inline std::complex<double> UProductFlipped(int& s,int& p,int kIndex);
        void setPmatrix(Eigen::MatrixXi& Pmatrix);
        int genNumbKMatrices(Eigen::MatrixXi& Pmatrix);
        double doublefactorial(int x);
        int g(int n,int m);
        void genKMatrix(Eigen::MatrixXi& Pmatrix,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col);
        Eigen::VectorXi genNonZeroCoordinates(int nonZeroEntries, Eigen::MatrixXi subTempKMatrix);
        void blockAssign(Eigen::MatrixXi tempKMatrix,int l,Eigen::MatrixXi& subTempKMatrix);
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        Eigen::MatrixXi genContributions(Eigen::VectorXi& impactedKet,Eigen::VectorXi& impactingKet,int& subPhotons, int& subModes);
        Eigen::MatrixXi validInputsAndOutputs(int& subPhotons, int& subModes,Eigen::VectorXi& in);
        int numbNonZeroEntries(Eigen::MatrixXi in);
        void genCijMatrix(Eigen::MatrixXi& Pmatrix,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col);
        void setLdimension(Eigen::MatrixXi& Pmatrix);
        Eigen::MatrixXcd genHermitian(Eigen::ArrayXd& a);
        Eigen::MatrixXcd genUnitary(Eigen::ArrayXd& a);
        Eigen::MatrixXcd matrixExp(Eigen::MatrixXcd X);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> saes;
        bool containedIn(int& ind,Eigen::ArrayXi& Arr);

};


inline int LOTransform::genKIndex(int& i,int& j,int& dim){
    return j + dim * i - (i*(i+1))/2;
}


inline std::complex<double> LOTransform::UProduct(int& s,int& p,int kIndex){
    std::complex<double> output(1,0);
    for(int q=0;q<kMatrix.at(kIndex)(p);q++){
        output *= pow(U(kMatrix.at(kIndex)(p+3*q+1),kMatrix.at(kIndex)(p+3*q+2)),kMatrix.at(kIndex)(p+3*q+3));
    }
    return output;
}


inline std::complex<double> LOTransform::UProductFlipped(int& s,int& p,int kIndex){
    std::complex<double> output(1,0);
    for(int q=0;q<kMatrix.at(kIndex)(p);q++){
        output *= pow(U(kMatrix.at(kIndex)(p+3*q+2),kMatrix.at(kIndex)(p+3*q+1)),kMatrix.at(kIndex)(p+3*q+3));
    }
    return output;
}

inline Eigen::MatrixXcd LOTransform::genUnitary(Eigen::ArrayXd& a){
    return matrixExp(genHermitian(a));
}

inline Eigen::MatrixXcd LOTransform::genHermitian(Eigen::ArrayXd& a){
    int Hsize = sqrt(a.size());
    Eigen::MatrixXcd m(Hsize,Hsize);
    int extractIndex=0;
    for(int i=0;i<Hsize;i++){
        m(i,i)=a(extractIndex);
        extractIndex++;
        for(int j=i;j<Hsize;j++){
            if(i!=j){
                //m(i,j)=a(extractIndex)+I*a(extractIndex+1);
                //m(j,i)=a(extractIndex)-I*a(extractIndex+1);
                m(i,j) = a(extractIndex) * exp(I*a(extractIndex+1));
                m(j,i) = a(extractIndex) * exp(-I*a(extractIndex+1));
                extractIndex++;
                extractIndex++;
            }
        }
    }


    return m;
}


inline Eigen::MatrixXcd LOTransform::matrixExp(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = exp(I*evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result=result+exp(I*evalues(j))*sylvester[j];
    }

    return result;
}

#endif // LOTRANSFORM_H_INCLUDED
