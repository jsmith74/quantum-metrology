#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>

/** ==== MAKE SURE THIS TOLERANCE IS ACCEPTABLE FOR GRID YOU USE ==== */

#define eps 1e-10

/** ================================================================= */

bool eqDoubs(double&d1,double& d2){

    if(d2 <= d1 + eps && d2 >= d1-eps) return true;

    return false;

}

void setMinandMax(std::string& preVarianceAddress,std::string& postVarianceAddress,double& max,double& min){

    std::ifstream infile(preVarianceAddress.c_str());

    while(!infile.eof()){

        double placeHold,throwaway;

        infile >> placeHold;

        infile >> throwaway;

        if(placeHold < min) min = placeHold;
        if(placeHold > max) max = placeHold;

    }

    infile.close();

    infile.open(postVarianceAddress.c_str());

    while(!infile.eof()){

        double placeHold,throwaway;

        infile >> placeHold;

        infile >> throwaway;

        if(placeHold < min) min = placeHold;
        if(placeHold > max) max = placeHold;

    }

    return;

}

double findMinInterval(std::string& fileAddress){

    double throwaway;

    double interval = 1e300;

    for(int i=1;true;i++){

        std::ifstream infile(fileAddress.c_str());

        double bnd1,bnd2;

        for(int j=0;j<i;j++){

            infile >> bnd1;
            infile >> throwaway;

        }

        if(infile.eof()){

            infile.close();

            break;

        }

        infile.close();

        infile.open(fileAddress.c_str());

        while(!infile.eof()){

            infile >> bnd2;
            infile >> throwaway;

            if(!eqDoubs(bnd1,bnd2)){

                double intervalHolder = std::abs(bnd1 - bnd2);

                if(intervalHolder < interval) interval = intervalHolder;

            }

        }

    }

    return interval;

}

void setIntervalSize(std::string& preVarianceAddress,std::string& postVarianceAddress,double& interval){

    interval = findMinInterval(preVarianceAddress);

    double interval2 = findMinInterval(postVarianceAddress);

    interval = std::max(interval,interval2);

    return;

}

bool fileSearch(double& delta,double& variance,std::string& fileAddress){

    std::ifstream infile(fileAddress.c_str());

    while(!infile.eof()){

        double deltaTest;

        infile >> deltaTest;

        infile >> variance;

        if(eqDoubs(delta,deltaTest)) return true;

    }

    infile.close();

    return false;

}

int main(int argc,char *argv[]){

    if(argc != 4){

        std::cout << "./VarianceCompare [Pre-variance File Address] [Post-variance File address] [% Improvement File Address]" << std::endl;

        return 1;

    }

    std::string preVarianceAddress,postVarianceAddress,improvementFileAddress;

    preVarianceAddress = argv[1];

    postVarianceAddress = argv[2];

    improvementFileAddress = argv[3];

    double min = 1e300;
    double max = -1e300;
    double interval;

    setMinandMax(preVarianceAddress,postVarianceAddress,max,min);

    setIntervalSize(preVarianceAddress,postVarianceAddress,interval);

    double variancePre,variancePost;

    std::ofstream outfile(improvementFileAddress.c_str());

    for(double i=min;i<=max;i+=interval){

        if(fileSearch(i,variancePre,preVarianceAddress) && fileSearch(i,variancePost,postVarianceAddress)){

            outfile << i << "\t" << (variancePre - variancePost) / variancePre << std::endl;

        }

    }

    outfile.close();

	return 0;

}
