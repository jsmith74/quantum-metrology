#include "BranchMeasStruct.h"

// TO DO: THINK OF A CLEVER WAY TO I/O A BRANCH STRUCTURE - WE WILL WANT IT TO BE AUTOMATIC
// TO DO: IF GENERATING BRANCH STRUCTURE BECOMES A BOTTLENECK, FIND A WAY TO SHARE LOOP OBJECTS
//        BETWEEN IDENTICAL STRUCTURES

BranchMeasStruct::BranchMeasStruct(){

}

void BranchMeasStruct::setKernalProbDistribution(){

    /** ======= ENTER THE INITIAL PROBABILITY DISTRIBUTION OF THE PHASE UNCERTAINTY ======================== */

        for(int i=0;i<numbGridPoints;i++){

            chainMeasurement.at(0).at(0).P_phi.at(i) = 1.0 / ( 2.0 * delta  );

        }

    /** ==================================================================================================== */

    chainMeasurement.at(0).at(0).printPDist();

    //Integration integrate;

    //integrate.setIntegral(delta,dP);

    //std::cout << std::setprecision(16)  <<  integrate.integrateArray(chainMeasurement.at(0).at(0).P_phi) << std::endl;

    return;

}


void BranchMeasStruct::setPhaseEstimators(){

    for(int i=0;i<numbTotalMeasBranches;i++){

        setBArray(i);
        printBArray();

        if(adaptive){

            for(int j=0;j<chainMeasurement[levels-1][i].numbBranches;j++){

                setMArrayAdaptive(j);
                printMArray();

                // UP TO HERE: CHECK THAT THE m AND b ARRAYS ARE GENERATED PROPERLY
                // THEN put code that uses m and b here ALSO ADD NON-ADAPTIVE VERSION

            }
            std::cout << std::endl;

        }

        else{

            for(int j=0;j<numbTotalMeasOutcomes;j++){



            }

            assert(1>2 && "TO DO: WRITE THIS");

        }

    }

    return;

}


inline void BranchMeasStruct::setBArray(int& i){

    b[levels-1] = i;

    for(int j=levels-1;j>0;j--){

        b[j-1] = chainMeasurement[j][ b[j] ].root;

    }

    return;

}


inline void BranchMeasStruct::setMArrayAdaptive(int& j){

    m[levels-1] = j;

    for(int i=levels-1;i>0;i--){

        m[i-1] = chainMeasurement[i][ b[i] ].rootMeas;

    }

    return;

}

void BranchMeasStruct::printMArray(){

    for(int i=0;i<m.size();i++){

        std::cout << m[i] << " ";

    }

    std::cout << std::endl;

    return;

}

void BranchMeasStruct::printBArray(){

    for(int i=0;i<b.size();i++){
        std::cout << b[i] << " ";
    }

    std::cout << std::endl;

    return;

}



void BranchMeasStruct::setPsiAndGamma(Eigen::VectorXd& position){

    int k=0;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement[i].size();j++){

            chainMeasurement[i][j].setPsi(position,k);

            chainMeasurement[i][j].updateGamma(position(k));
            k++;

        }

    }

    return;

}

int BranchMeasStruct::setFuncDimension(){

    int output = 0;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement.at(i).size();j++){

            output += 2*chainMeasurement.at(i).at(j).subHSDimension;

            output += 1;

        }

    }

    return output;

}

void BranchMeasStruct::printBranchStructure(){

    std::cout << "Levels: " << levels << std::endl << std::endl;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement.at(i).size();j++){


            std::cout << "level: " << chainMeasurement.at(i).at(j).level << std::endl;
            std::cout << "Index: " << j << std::endl;
            if(i>0) std::cout << "root: " << chainMeasurement.at(i).at(j).root << std::endl;
            if(i>0) std::cout << "rootMeas: " << chainMeasurement.at(i).at(j).rootMeas << std::endl;
            std::cout << "photons: " << chainMeasurement.at(i).at(j).extractPhotons() << std::endl;
            //if(i<levels-1){

                std::cout << "numbBranches: " << chainMeasurement.at(i).at(j).numbBranches << std::endl;

                std::cout << "branches:\n";
                for(int k=0;k<chainMeasurement.at(i).at(j).numbBranches;k++) std::cout << chainMeasurement.at(i).at(j).branches.at(k) << " ";

                std::cout << std::endl;

                chainMeasurement.at(i).at(j).printMAddress();

            //}

            std::cout << std::endl << std::endl;

        }

    }

    return;

}

void BranchMeasStruct::setKernel(){

    chainMeasurement.at(0).resize(1);

    chainMeasurement.at(0).at(0).initializeMZIObject(2,4,2);

    chainMeasurement.at(0).at(0).level = 0;

    chainMeasurement.at(0).at(0).delta = delta;

    chainMeasurement.at(0).at(0).dP = dP;

    chainMeasurement.at(0).at(0).P_phi.resize(numbGridPoints);

    return;

}

void BranchMeasStruct::setAdaptiveMeasurements(){

    for(int i=0;i<levels-1;i++){

        chainMeasurement.at(i+1).resize(0);

        int k=0;

        for(int j=0;j<chainMeasurement.at(i).size();j++){

            chainMeasurement.at(i+1).resize(chainMeasurement.at(i+1).size() + chainMeasurement.at(i).at(j).numbBranches);

            int l=0;

            while(k<chainMeasurement.at(i+1).size()){

                chainMeasurement.at(i+1).at(k).initializeMZIObject(2,4,2);    // THINK OF SOME WAY TO CHANGE THIS DYNAMICALLY

                chainMeasurement.at(i+1).at(k).root = j;

                chainMeasurement.at(i+1).at(k).rootMeas = l;

                chainMeasurement.at(i+1).at(k).level = i+1;

                chainMeasurement.at(i).at(j).branches.at(l) = k;

                chainMeasurement.at(i+1).at(k).delta = delta;

                chainMeasurement.at(i+1).at(k).dP = dP;

                //chainMeasurement.at(i+1).at(k).P_phi.resize(numbGridPoints);

                l++;

                k++;

            }

        }

    }

    return;

}

void BranchMeasStruct::setNonAdaptiveMeasurements(){

    for(int i=0;i<levels-1;i++){

        chainMeasurement.at(i+1).resize(1);

        chainMeasurement.at(i+1).at(0).initializeMZIObject(2,4,2);

        chainMeasurement.at(i+1).at(0).root = 0;

        chainMeasurement.at(i+1).at(0).level = i+1;

        chainMeasurement.at(i+1).at(0).delta = delta;

        chainMeasurement.at(i+1).at(0).dP = dP;

        //chainMeasurement.at(i+1).at(0).P_phi.resize(numbGridPoints);

        for(int j=0;j<chainMeasurement.at(i).at(0).numbBranches;j++) chainMeasurement.at(i).at(0).branches.at(j) = 0;

    }

    return;

}

void BranchMeasStruct::setNumbTotalMeasOutcomesAdaptive(){

    numbTotalMeasOutcomes = 0;

    for(int i=0;i<numbTotalMeasBranches;i++){

        numbTotalMeasOutcomes += chainMeasurement.at(levels-1).at(i).numbBranches;

    }

    return;

}

void BranchMeasStruct::setNumbTotalMeasOutcomesNonAdaptive(){

    numbTotalMeasOutcomes = 1;

    for(int i=0;i<levels;i++){

        numbTotalMeasOutcomes *= chainMeasurement.at(i).at(0).numbBranches;

    }

    return;

}

void BranchMeasStruct::setMeasChain(bool Adaptive,int numbMeas,bool Import,int gridSize,double Delta){

    adaptive = Adaptive;
    levels = numbMeas;
    import = Import;

    delta = Delta;

    dP = (2.0 * delta) / (1.0 * gridSize);

    numbGridPoints = gridSize + 1;

    assert(numbGridPoints % 2 == 1 && "ERROR: Simpson's rule requires an even number of intervals (odd number of grid points).");

    chainMeasurement.resize(levels);

    setKernel();

    if(adaptive) setAdaptiveMeasurements();

    if(!adaptive) setNonAdaptiveMeasurements();

    numbTotalMeasBranches = chainMeasurement.at(levels-1).size();

    if(adaptive) setNumbTotalMeasOutcomesAdaptive();
    if(!adaptive) setNumbTotalMeasOutcomesNonAdaptive();

    phaseEstimators.resize(numbTotalMeasOutcomes);

    m.resize(levels);

    b.resize(levels);

    return;

}
