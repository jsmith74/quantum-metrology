#include "PUA.h"


//#define PRINT_BASIS_CHECK


/** ====== Restrict Unitarity of Linear Operators (also change in MeritFunction.cpp and LOTransform.cpp  ============= */

#define RESTRICT_TO_UNITARY

//#define ALLOW_ARBITRARY_VACUUM_MODES

/** ================================================================================================================== */


void PUA::setQuantumOperator(Eigen::VectorXd& position){

//todo: speed this up; eliminate requirement that we copy these vectors
    Eigen::ArrayXd LOOPPosition = position.segment(0,LOOPPositionSize);
    Eigen::ArrayXd AAugPosition = position.segment(LOOPPositionSize,AAugPositionSize);

    LOOP.setUnitaryMatrix(LOOPPosition);
    AAug.setAugmentMatrix(AAugPosition);

    //AAug.printAugmentMatrix();

    for(int i=0;i<nonZeroRows;i++){

        for(int j=0;j<nonZeroCols;j++){

            OMEGAU(i,j) = LOOP.omegaUij(Row(i),Col(j));

        }

    }


    TotalOp = OMEGAU * AAug.AugmentMatrix;

    return;

}

PUA::PUA(){

}

PUA::PUA(int N,int M,int Na,int Ma,int Mm,int outcome){

    AAug.setAncillaAugment(N,M,Na,Ma);

    #ifdef PRINT_BASIS_CHECK

        AAug.printAugmentMatrix();

    #endif // PRINT_BASIS_CHECK

    Eigen::MatrixXi initialBasisVector = generateBasisVector(N,M,1);
    Eigen::MatrixXi midBasisVector = generateBasisVector(N+Na,M+Ma,1);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Initial Basis Vec:\n" << initialBasisVector << std::endl << std::endl;
        std::cout << "mid basis vec: \n" << midBasisVector << std::endl << std::endl;

    #endif // PRINT_BASIS_CHECK

    measOutcome.resize(Mm);

    setMeasOutcome(N+Na,M+Ma,Mm,outcome);

    std::ofstream outfile("Measurement_Outcome.dat");
    outfile << measOutcome << std::endl;
    outfile.close();

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Meas Outcome Vec:\n" << measOutcome << std::endl << std::endl;

    #endif

    int NLeft = N + Na - measOutcome.sum();
    int MLeft = M + Ma - Mm;

    #ifdef PRINT_BASIS_CHECK

        std::cout << "NLeft: " << NLeft << std::endl;
        std::cout << "MLeft: " << MLeft << std::endl;

    #endif

    Eigen::MatrixXi finalBasisVector = generateBasisVector(NLeft,MLeft,1);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Final Basis Vec:\n" << finalBasisVector << std::endl << std::endl;

    #endif

    nonZeroCols = g(N,M) * g(Na,Ma);
    nonZeroRows = g(NLeft,MLeft);

    if(Na==0 && Ma==0) nonZeroCols = g(N,M);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "NonZeroCols: " << nonZeroCols << std::endl;
        std::cout << "NonZeroRows: " << nonZeroRows << std::endl;

    #endif

    Col.resize(nonZeroCols);
    Row.resize(nonZeroRows);

    TotalOp.resize(nonZeroRows,g(N,M));

    setCol(N,M,Na,Ma,midBasisVector,initialBasisVector);

    setRow(midBasisVector,finalBasisVector);

    LOOP.setLOTransform(N+Na,M+Ma,Row,Col);

    #ifdef PRINT_BASIS_CHECK

        for(int i=0;i<Col.size();i++) std::cout << midBasisVector.row(Col(i)) << std::endl;
        std::cout << std::endl;
        for(int i=0;i<Row.size();i++) std::cout << midBasisVector.row(Row(i)) << std::endl;

    #endif // PRINT_BASIS_CHECK

    #ifdef RESTRICT_TO_UNITARY

        LOOPPositionSize = (M + Ma) * (M + Ma);

    #endif // RESTRICT_TO_UNITARY

    #ifdef ALLOW_ARBITRARY_VACUUM_MODES

        LOOPPositionSize = 2 * (M + Ma) * (M + Ma);

    #endif // ALLOW_ARBITRARY_VACUUM_MODES

    AAugPositionSize = 2*g(Na,Ma);

    OMEGAU.resize(nonZeroRows,nonZeroCols);

}


void PUA::setPUA(int N,int M,int Na,int Ma,int Mm,int outcome){

    AAug.setAncillaAugment(N,M,Na,Ma);

    #ifdef PRINT_BASIS_CHECK

        AAug.printAugmentMatrix();

    #endif // PRINT_BASIS_CHECK

    Eigen::MatrixXi initialBasisVector = generateBasisVector(N,M,1);
    Eigen::MatrixXi midBasisVector = generateBasisVector(N+Na,M+Ma,1);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Initial Basis Vec:\n" << initialBasisVector << std::endl << std::endl;
        std::cout << "mid basis vec: \n" << midBasisVector << std::endl << std::endl;

    #endif

    measOutcome.resize(Mm);

    setMeasOutcome(N+Na,M+Ma,Mm,outcome);

    std::ofstream outfile("Measurement_Outcome.dat");
    outfile << measOutcome << std::endl;
    outfile.close();

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Meas Outcome Vec:\n" << measOutcome << std::endl << std::endl;

    #endif

    int NLeft = N + Na - measOutcome.sum();
    int MLeft = M + Ma - Mm;

    #ifdef PRINT_BASIS_CHECK

        std::cout << "NLeft: " << NLeft << std::endl;
        std::cout << "MLeft: " << MLeft << std::endl;

    #endif

    Eigen::MatrixXi finalBasisVector = generateBasisVector(NLeft,MLeft,1);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "Final Basis Vec:\n" << finalBasisVector << std::endl << std::endl;

    #endif // PRINT_BASIS_CHECK

    nonZeroCols = g(N,M) * g(Na,Ma);
    nonZeroRows = g(NLeft,MLeft);

    if(Na==0 && Ma==0) nonZeroCols = g(N,M);

    #ifdef PRINT_BASIS_CHECK

        std::cout << "NonZeroCols: " << nonZeroCols << std::endl;
        std::cout << "NonZeroRows: " << nonZeroRows << std::endl;

    #endif

    Col.resize(nonZeroCols);
    Row.resize(nonZeroRows);

    TotalOp.resize(nonZeroRows,g(N,M));

    setCol(N,M,Na,Ma,midBasisVector,initialBasisVector);

    setRow(midBasisVector,finalBasisVector);

    LOOP.setLOTransform(N+Na,M+Ma,Row,Col);

    #ifdef PRINT_BASIS_CHECK

        for(int i=0;i<Col.size();i++) std::cout << midBasisVector.row(Col(i)) << std::endl;
        std::cout << std::endl;
        for(int i=0;i<Row.size();i++) std::cout << midBasisVector.row(Row(i)) << std::endl;

    #endif // PRINT_BASIS_CHECK

    #ifdef RESTRICT_TO_UNITARY

        LOOPPositionSize = (M + Ma) * (M + Ma);

    #endif // RESTRICT_TO_UNITARY

    #ifdef ALLOW_ARBITRARY_VACUUM_MODES

        LOOPPositionSize = 2 * (M + Ma) * (M + Ma);

    #endif // ALLOW_ARBITRARY_VACUUM_MODES

    AAugPositionSize = 2*g(Na,Ma);

    OMEGAU.resize(nonZeroRows,nonZeroCols);

}


void PUA::setRow(Eigen::MatrixXi& midBasisVector,Eigen::MatrixXi& finalBasisVector){

    int k=0;

    for(int i=0;i<finalBasisVector.rows();i++){

        Eigen::ArrayXi testArray1(midBasisVector.cols());

        testArray1.segment(0,measOutcome.size()) = measOutcome.transpose();

        testArray1.segment(measOutcome.size(),finalBasisVector.cols()) = finalBasisVector.row(i).transpose();

        for(int j=0;j<midBasisVector.rows();j++){

            Eigen::ArrayXi testArray2(midBasisVector.cols());
            testArray2 = midBasisVector.row(j).transpose();

            if(isEqual(testArray1,testArray2)){

                Row(k) = j;
                k++;
                break;

            }

        }


    }

    return;

}

void PUA::setCol(int N,int M,int Na,int Ma,Eigen::MatrixXi& midBasisVector,Eigen::MatrixXi& initialBasisVector){

    if(Na==0 && Ma==0){
        for(int i=0;i<Col.size();i++) Col(i) = i;
        return;
    }

    Eigen::MatrixXi ancillaBasis = generateBasisVector(Na,Ma,1);

    int k=0;

    for(int i=0;i<ancillaBasis.rows();i++){

        for(int j=0;j<initialBasisVector.rows();j++){

            Eigen::ArrayXi testArray1(M+Ma);
            testArray1.segment(0,Ma) = ancillaBasis.row(i).transpose();
            testArray1.segment(Ma,M) = initialBasisVector.row(j).transpose();

            for(int l=0;l<midBasisVector.rows();l++){

                Eigen::ArrayXi testArray2(M+Ma);
                testArray2 = midBasisVector.row(l).transpose();

                if(isEqual(testArray1,testArray2)){

                    Col(k) = l;
                    k++;
                    break;

                }

            }

        }

    }


    return;

}

bool PUA::isEqual(Eigen::ArrayXi a1,Eigen::ArrayXi a2){

    assert(a1.size() == a2.size() && "Error: arrays of different sizes are being compared.");

    for(int i=0;i<a1.size();i++){

        if(a1(i) != a2(i)) return false;

    }

    return true;

}

void PUA::setMeasOutcome(int N,int M,int Mm,int outcome){

    if(Mm==0) return;

    int totalRows = 0;

    for(int Nm=0;Nm<=N;Nm++){

        Eigen::MatrixXi basisVector = generateBasisVector(Nm,Mm,1);

        totalRows += basisVector.rows();

        if(outcome < totalRows){

            totalRows -= basisVector.rows();

            measOutcome = basisVector.row(outcome - totalRows);

            break;

        }

    }

    return;

}

Eigen::MatrixXi PUA::generateSubBasisVector(int subPhotons, int subModes){

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



Eigen::MatrixXi PUA::generateBasisVector(int subPhotons,int subModes, int subMeasureModes){
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


int PUA::g(int n,int m){
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


double PUA::doublefactorial(int x){
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
