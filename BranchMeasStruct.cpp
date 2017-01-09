#include "BranchMeasStruct.h"


BranchMeasStruct::BranchMeasStruct(){

}

void BranchMeasStruct::printBranchStructure(){

    std::cout << "Levels: " << levels << std::endl << std::endl;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement.at(i).size();j++){

            std::cout << i << " " << j << std::endl;
            std::cout << "photons: " << chainMeasurement.at(i).at(j).extractPhotons() << std::endl;
            std::cout << "level: " << chainMeasurement.at(i).at(j).level << std::endl;

            if(i>0) std::cout << "root: " << chainMeasurement.at(i).at(j).root << std::endl;

            if(i<levels-1){

                std::cout << "numbBranches: " << chainMeasurement.at(i).at(j).numbBranches << std::endl;

                std::cout << "branches:\n";
                for(int k=0;k<chainMeasurement.at(i).at(j).numbBranches;k++) std::cout << chainMeasurement.at(i).at(j).branches.at(k) << " ";

                std::cout << std::endl;

            }

            std::cout << std::endl;

        }

    }

    return;

}

void BranchMeasStruct::setKernel(){

    chainMeasurement.at(0).resize(1);

    chainMeasurement.at(0).at(0).initializeMZIObject(2,4,2);

    chainMeasurement.at(0).at(0).level = 0;

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

                chainMeasurement.at(i+1).at(k).initializeMZIObject(l%4+1,4,2);    // THINK OF SOME WAY TO CHANGE THIS DYNAMICALLY

                chainMeasurement.at(i+1).at(k).root = j;

                chainMeasurement.at(i+1).at(k).level = i+1;

                chainMeasurement.at(i).at(j).branches.at(l) = k;

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

        for(int j=0;j<chainMeasurement.at(i).at(0).numbBranches;j++) chainMeasurement.at(i).at(0).branches.at(j) = 0;

    }

    return;

}

void BranchMeasStruct::setMeasChain(bool Adaptive,int numbMeas,bool Import){

    adaptive = Adaptive;
    levels = numbMeas;
    import = Import;

    chainMeasurement.resize(levels);

    setKernel();

    if(adaptive) setAdaptiveMeasurements();

    if(!adaptive) setNonAdaptiveMeasurements();

    return;

}
