#include "LOTransform.h"


std::complex<double> LOTransform::omegaUij(int& i,int& j){

    if(i<=j){
        int kIndex = genKIndex(i,j,LDimension(1));
        std::complex<double> output(0,0);

        int p=1;
        for(int q=0;q<kMatrix[kIndex](0);q++){
            output += CijMatrix[kIndex](q) * UProduct(q,p,kIndex);
            p = p + 3*kMatrix[kIndex](p) + 1;
        }
        return output;
    }

    else{
        int kIndex = genKIndex(j,i,LDimension(1));
        std::complex<double> output(0,0);
        int p=1;
        for(int q=0;q<kMatrix[kIndex](0);q++){
            output += CijMatrix[kIndex](q) * UProductFlipped(q,p,kIndex);
            p = p + 3*kMatrix[kIndex](p) + 1;
        }
        return output;
    }

}

void LOTransform::setUnitaryMatrixDirect(Eigen::MatrixXcd& Uin){

    U = Uin;

    return;

}

void LOTransform::setUnitaryMatrix(Eigen::ArrayXd& position){

    #ifdef RESTRICT_TO_UNITARY

        U = genUnitary(position);

    #endif // RESTRICT_TO_UNITARY

    #ifdef ALLOW_ARBITRARY_VACUUM_MODES

        int k=0;

        for(int i=0;i<modes;i++){

            for(int j=0;j<modes;j++){

                U(i,j) = position(k);
                k++;

                U(i,j) += position(k) * I;
                k++;

            }

        }

        saes.compute(U.conjugate().transpose() * U);

        U *= 1.0/sqrt(saes.eigenvalues().maxCoeff());

    #endif // ALLOW_ARBITRARY_VACUUM_MODES

    return;
}

LOTransform::LOTransform(){


}

LOTransform::LOTransform(int N,int M,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col){

    photons = N;
    modes = M;

    std::complex<double> IGen(0.0,1.0);
    I = IGen;

    Eigen::MatrixXi Pmatrix(3,3);

    setPmatrix(Pmatrix);

    const int numbKMatrices = genNumbKMatrices(Pmatrix);

    kMatrix.resize(numbKMatrices);
    genKMatrix(Pmatrix,Row,Col);

    CijMatrix.resize(numbKMatrices);
    genCijMatrix(Pmatrix,Row,Col);

    setLdimension(Pmatrix);

    U.resize(modes,modes);

}


void LOTransform::setLOTransform(int N,int M,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col){

    photons = N;
    modes = M;

    std::complex<double> IGen(0.0,1.0);
    I = IGen;

    Eigen::MatrixXi Pmatrix(3,3);

    setPmatrix(Pmatrix);

    const int numbKMatrices = genNumbKMatrices(Pmatrix);

    kMatrix.resize(numbKMatrices);
    genKMatrix(Pmatrix,Row,Col);

    CijMatrix.resize(numbKMatrices);
    genCijMatrix(Pmatrix,Row,Col);

    setLdimension(Pmatrix);

    U.resize(modes,modes);

}


void LOTransform::setLdimension(Eigen::MatrixXi& Pmatrix){

    int inPhotons = Pmatrix(2,0) - Pmatrix(1,0);
    int inModes = Pmatrix(0,1) - Pmatrix(0,2);
    int outPhotons = Pmatrix(2,0);
    int outModes = Pmatrix(2,1);
    LDimension << g(inPhotons,inModes),g(outPhotons,outModes);

    return;

}


void LOTransform::genCijMatrix(Eigen::MatrixXi& Pmatrix,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col){

    int dim = g(Pmatrix(2,0),Pmatrix(2,1));

    for(int i=0;i<dim;i++){
        for(int j=i;j<dim;j++){

            if((containedIn(i,Row) &&  containedIn(j,Col)) || (containedIn(i,Col) &&  containedIn(j,Row))){
                CijMatrix[genKIndex(i,j,dim)] = Eigen::ArrayXd::Zero(kMatrix[genKIndex(i,j,dim)](0));
            }

        }
    }

    Eigen::MatrixXi basisVector = generateBasisVector(Pmatrix(2,0),Pmatrix(2,1),Pmatrix(2,2));

    for(int i=0;i<dim;i++){
        for(int j=i;j<dim;j++){

            if((containedIn(i,Row) &&  containedIn(j,Col)) || (containedIn(i,Col) &&  containedIn(j,Row))){

                for(int l=0;l<kMatrix[genKIndex(i,j,dim)](0);l++){
                    Eigen::VectorXi impactingKet = basisVector.row(j);
                    Eigen::VectorXi impactedKet = basisVector.row(i);
                    CijMatrix[genKIndex(i,j,dim)](l) = 1;
                    for(int q=0;q<Pmatrix(2,1);q++){
                        CijMatrix[genKIndex(i,j,dim)](l) = CijMatrix[genKIndex(i,j,dim)](l)*doublefactorial(impactingKet(q))*doublefactorial(impactedKet(q));
                    }
                    CijMatrix[genKIndex(i,j,dim)](l) = sqrt(CijMatrix[genKIndex(i,j,dim)](l));
                }
            }

        }
    }

    for(int i=0;i<dim;i++){
        for(int j=i;j<dim;j++){
            if((containedIn(i,Row) &&  containedIn(j,Col)) || (containedIn(i,Col) &&  containedIn(j,Row))){
                Eigen::ArrayXi tempKMatrix = kMatrix[genKIndex(i,j,dim)];
                int p = 1;
                for(int l=0;l<tempKMatrix(0);l++){
                    for(int q=0;q<tempKMatrix(p);q++){
                        CijMatrix[genKIndex(i,j,dim)](l) = CijMatrix[genKIndex(i,j,dim)](l) * 1.0/doublefactorial(tempKMatrix(p+3*(q+1)));
                    }
                    p = p + 3*tempKMatrix(p) + 1;
                }
            }
        }
    }
    return;
}



Eigen::MatrixXi LOTransform::generateSubBasisVector(int subPhotons, int subModes){
    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    Eigen::MatrixXi nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );
    return nv;
}


int LOTransform::numbNonZeroEntries(Eigen::MatrixXi in){

    int output = 0;
    for(int i=0;i<in.rows();i++){
        for(int j=0;j<in.cols();j++){
            if(in(i,j)!=0){
                output++;
            }
        }
    }

    return output;

}

Eigen::MatrixXi LOTransform::genContributions(Eigen::VectorXi& impactedKet,Eigen::VectorXi& impactingKet,int& subPhotons, int& subModes){

    Eigen::MatrixXi tempValidInputsandOutputs = validInputsAndOutputs(subPhotons,subModes,impactingKet);
    int numbSolutions = 0;
    int numbCombinations = 1;
    for(int i=0;i<subModes;i++){
        numbCombinations = numbCombinations*tempValidInputsandOutputs(0,i);
    }
    Eigen::VectorXi numbStorage(subModes);
    numbStorage = tempValidInputsandOutputs.row(0);

    Eigen::VectorXi numbIteration = Eigen::VectorXi::Zero(subModes);
    Eigen::MatrixXi validSolutions(0,subModes);
    for(int i=0;i<numbCombinations;i++){
        Eigen::MatrixXi tempCombinationStorage(subModes,subModes);
        for(int j=0;j<subModes;j++){
            int fixIndex =0;
            for(int k=0;k<j;k++){
                fixIndex = fixIndex + numbStorage(k);
            }
            tempCombinationStorage.row(j) = tempValidInputsandOutputs.row(1+numbIteration(j)+fixIndex);
        }

        //append it to solutions if it checks out

        numbIteration(0)++;
        for(int j=0;j<subModes;j++){
            bool increTest = false;
            if(numbIteration(j) > 0){
                increTest = true;
            }
            numbIteration(j) = numbIteration(j) % numbStorage(j);
            if(numbIteration(j)==0 && j < subModes-1 && increTest){
                numbIteration(j+1)++;
            }
        }

        for(int j=0;j<subModes;j++){
            if((tempCombinationStorage.col(j).sum())!=impactedKet(j)){
                break;
            }
            if(j==subModes-1){
                numbSolutions++;
                validSolutions.conservativeResize(subModes*numbSolutions,subModes);

                for(int k=0;k<subModes;k++){
                    for(int l=0;l<subModes;l++){
                        validSolutions(k+subModes*(numbSolutions-1),l) = tempCombinationStorage(k,l);
                    }
                }
            }
        }

    }

    return validSolutions;

}


Eigen::MatrixXi LOTransform::validInputsAndOutputs(int& subPhotons, int& subModes,Eigen::VectorXi& in){
    Eigen::MatrixXi outMatrix = Eigen::MatrixXi::Zero(1,subModes);
    int numbSolutions = 1;
    for(int j=0; j<subModes;j++){
        Eigen::VectorXi holdVector = Eigen::VectorXi::Zero(subModes);
        int numbComb = 1;
        for(int i =0;i<subModes;i++){
            numbComb=numbComb*(in(j)+1);
        }
        for(int i = 0;i<numbComb;i++){
            if(holdVector.sum()==in(j)){

                outMatrix.conservativeResize(numbSolutions+1,subModes);

                outMatrix.row(numbSolutions) = holdVector;
                outMatrix(0,j)++;
                numbSolutions++;
            }
            holdVector(0)++;
            for(int k=0;k<subModes;k++){
                if(holdVector(k)>in(j)){
                    holdVector(k) = 0;
                    if(k<subModes-1){
                        holdVector(k+1)++;
                    }
                }
            }
        }

    }

    return outMatrix;

}


Eigen::MatrixXi LOTransform::generateBasisVector(int subPhotons,int subModes, int subMeasureModes){
    Eigen::MatrixXi output(0,subModes);
    for(int i=subPhotons;i>=0;i--){

        Eigen::MatrixXi measModes = generateSubBasisVector(i,subMeasureModes);
        if(subMeasureModes==subModes){
            return measModes;
        }

        Eigen::MatrixXi othaModes = generateSubBasisVector(subPhotons-i,subModes-subMeasureModes);

        int numbRows = measModes.rows()*othaModes.rows();
        int numbMeasRows = measModes.rows();
        int numbOthaRows = othaModes.rows();
        int outputrows = output.rows();
        output.conservativeResize(outputrows+numbRows,subModes);

        for(int k=0;k<numbMeasRows;k++){
            for(int j=0;j<numbOthaRows;j++){
                for(int l=0;l<subMeasureModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = measModes(k,l);
                }
                for(int l=subMeasureModes;l<subModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = othaModes(j,l-subMeasureModes);
                }
            }
        }
    }
    return output;
}


void LOTransform::blockAssign(Eigen::MatrixXi tempKMatrix,int l,Eigen::MatrixXi& subTempKMatrix){
    int outModes = subTempKMatrix.rows();
    for(int i=0;i<outModes;i++){
        for(int j=0;j<outModes;j++){
            subTempKMatrix(i,j) = tempKMatrix(i+outModes*l,j);
        }
    }
    return;
}


Eigen::VectorXi LOTransform::genNonZeroCoordinates(int nonZeroEntries, Eigen::MatrixXi subTempKMatrix){

    Eigen::VectorXi outPut(3*nonZeroEntries);
    int k=0;
    for(int i=0;i<subTempKMatrix.rows();i++){
        for(int j=0;j<subTempKMatrix.cols();j++){

            if(subTempKMatrix(i,j) != 0){
                outPut(3*k) = i;
                outPut(3*k+1) = j;
                outPut(3*k+2) = subTempKMatrix(i,j);
                k++;
            }

        }
    }

    return outPut;

}

bool LOTransform::containedIn(int& ind,Eigen::ArrayXi& Arr){

    for(int i=0;i<Arr.size();i++){

        if(ind == Arr(i)) return true;

    }

    return false;

}

void LOTransform::genKMatrix(Eigen::MatrixXi& Pmatrix,Eigen::ArrayXi& Row,Eigen::ArrayXi& Col){

    int outPhotons = Pmatrix(2,0);
    int outModes = Pmatrix(2,1);
    int outMeasModes = Pmatrix(2,2);
    Eigen::MatrixXi basisVector = generateBasisVector(outPhotons,outModes,outMeasModes);
    int numbKArrays = (basisVector.rows()*basisVector.rows() + basisVector.rows())/2;
    for(int i=0;i<numbKArrays;i++){
        kMatrix.at(i) = Eigen::ArrayXi::Zero(2);
    }
    int k=0;

    for(int i=0;i<basisVector.rows();i++){
        for(int j=i;j<basisVector.rows();j++){

            if((containedIn(i,Row) &&  containedIn(j,Col)) || (containedIn(i,Col) &&  containedIn(j,Row))){


                Eigen::VectorXi impactingKet = basisVector.row(j);
                Eigen::VectorXi impactedKet = basisVector.row(i);
                Eigen::MatrixXi tempKmatrix = genContributions(impactedKet,impactingKet,outPhotons,outModes);
                int numbKMatrices = tempKmatrix.rows()/outModes;
                kMatrix.at(k)(0) = numbKMatrices;
                int p=1;
                for(int l=0;l<numbKMatrices;l++){

                    Eigen::MatrixXi subTempKMatrix(outModes,outModes);
                    blockAssign(tempKmatrix,l,subTempKMatrix);

                    int nonZeroEntries = numbNonZeroEntries(subTempKMatrix);
                    kMatrix.at(k)(p) = nonZeroEntries;
                    p++;

                    Eigen::VectorXi nonZeroCoordinates = genNonZeroCoordinates(nonZeroEntries,subTempKMatrix);

                    if(p + 3*nonZeroEntries >= kMatrix.at(k).size()){
                        int outputCols = kMatrix.at(k).size();
                        int difference = p+3*nonZeroEntries-outputCols;
                        kMatrix.at(k).conservativeResize(outputCols+difference+1);
                    }

                    for(int q=0;q<3*nonZeroEntries;q++){
                        kMatrix.at(k)(p+q) = nonZeroCoordinates(q);
                    }

                    p = p + 3*nonZeroEntries;

                }

            }

            k++;
        }
    }

    return;

}



int LOTransform::g(int n,int m){
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

double LOTransform::doublefactorial(int x){
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

int LOTransform::genNumbKMatrices(Eigen::MatrixXi& Pmatrix){
    int outPhotons = Pmatrix(2,0);
    int outModes = Pmatrix(2,1);
    return (g(outPhotons,outModes) * g(outPhotons,outModes) + g(outPhotons,outModes))/2;
}


void LOTransform::setPmatrix(Eigen::MatrixXi& Pmatrix){

    Pmatrix(0,0) = photons;
    Pmatrix(0,1) = modes;
    Pmatrix(0,2) = 0;
    Pmatrix(1,0) = 0;
    Pmatrix(1,1) = 0;
    Pmatrix(1,2) = 0;
    Pmatrix(2,0) = photons;
    Pmatrix(2,1) = modes;
    Pmatrix(2,2) = 1;
    return;

}


