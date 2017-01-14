#include "BranchMeasStruct.h"

// TO DO: THINK OF A CLEVER WAY TO I/O A BRANCH STRUCTURE - WE WILL WANT IT TO BE AUTOMATIC
// TO DO: IF GENERATING BRANCH STRUCTURE BECOMES A BOTTLENECK, FIND A WAY TO SHARE LOOP OBJECTS
//        BETWEEN IDENTICAL STRUCTURES

BranchMeasStruct::BranchMeasStruct(){

}

void BranchMeasStruct::setKernalProbDistribution(){

    /** ======= ENTER THE INITIAL PROBABILITY DISTRIBUTION OF THE PHASE UNCERTAINTY ======================== */

        for(int i=0;i<numbGridPoints;i++){

            chainMeasurement.at(0).at(0).P.at(i) = 1.0 / ( 2.0 * delta  );

        }

    /** ==================================================================================================== */

    //chainMeasurement.at(0).at(0).printPDist();

    return;

}

void BranchMeasStruct::setPhaseEstimators(){

    for(int i=0;i<levels;i++) m[i] = 0;

    // UP TO HERE, WRITE CODE THAT CONSTRUCTS PHASE ESTIMATORS GOING TO HAVE TO USE SIMPSONS RULE ETC
    // WRITE CODE TO DO SIMPSON'S RULE FIRST

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

    chainMeasurement.at(0).at(0).P.resize(numbGridPoints);

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

                chainMeasurement.at(i+1).at(k).level = i+1;

                chainMeasurement.at(i).at(j).branches.at(l) = k;

                chainMeasurement.at(i+1).at(k).delta = delta;

                chainMeasurement.at(i+1).at(k).dP = dP;

                chainMeasurement.at(i+1).at(k).P.resize(numbGridPoints);

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

        chainMeasurement.at(i+1).at(0).P.resize(numbGridPoints);

        for(int j=0;j<chainMeasurement.at(i).at(0).numbBranches;j++) chainMeasurement.at(i).at(0).branches.at(j) = 0;

    }

    return;

}

void BranchMeasStruct::setNumbTotalMeasOutcomesAdaptive(){

    numbTotalMeasOutcomes = 0;

    for(int i=0;i<chainMeasurement.at(levels-1).size();i++){

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

    numbGridPoints = gridSize;

    chainMeasurement.resize(levels);

    setKernel();

    if(adaptive) setAdaptiveMeasurements();

    if(!adaptive) setNonAdaptiveMeasurements();

    numbTotalMeasBranches = chainMeasurement.at(levels-1).size();

    if(adaptive) setNumbTotalMeasOutcomesAdaptive();
    if(!adaptive) setNumbTotalMeasOutcomesNonAdaptive();

    phaseEstimators.resize(numbTotalMeasOutcomes);

    m.resize(levels);

    return;

}
