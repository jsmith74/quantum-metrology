#include "BranchMeasStruct.h"

#define PI 3.141592653589793

#define PHOTONS 3


BranchMeasStruct::BranchMeasStruct(){



}


double BranchMeasStruct::generalVariance(){

    return integrate.generalVariance(chainMeasurement);

}


void BranchMeasStruct::setKernalProbDistribution(){

    /** ======= ENTER THE INITIAL PROBABILITY DISTRIBUTION OF THE PHASE UNCERTAINTY ======================== */

        for(int i=0;i<numbGridPoints;i++){

            chainMeasurement.at(0).at(0).P_phi.at(i) = 1.0 / ( 2.0 * delta  );

        }

    /** ==================================================================================================== */

    return;

}


void BranchMeasStruct::setMeasChain(bool Adaptive,int numbMeas,bool Import,int gridSize,double Delta){

    adaptive = Adaptive;
    levels = numbMeas;
    import = Import;

    delta = Delta;

    dP = (2.0 * delta) / (1.0 * gridSize);

    numbGridPoints = gridSize + 1;

    assert(numbGridPoints % 2 == 1 && "Simpson's rule requires an even number of intervals (odd number of grid points).");

    chainMeasurement.resize(levels);

    setKernel();

    if(adaptive) setAdaptiveMeasurements();

    if(!adaptive) setNonAdaptiveMeasurements();

    numbTotalMeasBranches = chainMeasurement.at(levels-1).size();

    if(adaptive) setNumbTotalMeasOutcomesAdaptive();
    if(!adaptive) setNumbTotalMeasOutcomesNonAdaptive();

    integrate.setIntegral(delta,dP,numbGridPoints,levels,adaptive,numbTotalMeasBranches,numbTotalMeasOutcomes);

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

void BranchMeasStruct::printStateAmps(Eigen::VectorXd& position,std::ofstream& outfile){

    int k=0;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement[i].size();j++){

            chainMeasurement[i][j].printStateAmps(position,k,outfile);

            k++;

        }

    }

    return;

}

void BranchMeasStruct::printGammaAmps(Eigen::VectorXd& position,std::ofstream& outfile){

    int k=0;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement[i].size();j++){

            k += 2 * chainMeasurement[i][j].subHSDimension;

            outfile << position(k) << "\t";

            k++;

        }

    }

    return;

}

void BranchMeasStruct::printPsiAndGamma(Eigen::VectorXd& position,std::ofstream& outfile){

    int k=0;

    for(int i=0;i<levels;i++){

        for(int j=0;j<chainMeasurement[i].size();j++){

            outfile << "Level: " << i << "\n";
            if(i>0) outfile << "Root: " << chainMeasurement[i][j].root << "\n";
            outfile << "Photons: " << chainMeasurement[i][j].extractPhotons() << std::endl << std::endl;

            chainMeasurement[i][j].printPsi(position,k,outfile);

            outfile << std::endl;

            outfile << "Gamma: " << std::setprecision(16) << position(k) << std::endl << std::endl << std::endl;

            k++;

        }

    }

    return;

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

            std::cout << "numbBranches: " << chainMeasurement.at(i).at(j).numbBranches << std::endl;

            std::cout << std::endl;

            chainMeasurement.at(i).at(j).printMAddress();

            std::cout << std::endl << std::endl;

        }

    }

    return;

}


void BranchMeasStruct::setKernel(){

    chainMeasurement.at(0).resize(1);

    chainMeasurement.at(0).at(0).initializeMZIObject(PHOTONS,2,2);

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

                chainMeasurement.at(i+1).at(k).initializeMZIObject(PHOTONS,2,2);

                chainMeasurement.at(i+1).at(k).root = j;

                chainMeasurement.at(i+1).at(k).rootMeas = l;

                chainMeasurement.at(i+1).at(k).level = i+1;

                chainMeasurement.at(i+1).at(k).delta = delta;

                chainMeasurement.at(i+1).at(k).dP = dP;

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

        chainMeasurement.at(i+1).at(0).initializeMZIObject(PHOTONS,2,2);

        chainMeasurement.at(i+1).at(0).root = 0;

        chainMeasurement.at(i+1).at(0).level = i+1;

        chainMeasurement.at(i+1).at(0).delta = delta;

        chainMeasurement.at(i+1).at(0).dP = dP;

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

