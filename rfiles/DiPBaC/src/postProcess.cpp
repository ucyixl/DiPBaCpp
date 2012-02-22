#include <vector>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "include/postProcess.h"
using std::string;
using std::ifstream;
using std::vector;

SEXP calcDisSimMat(SEXP fileName, SEXP nSweeps, SEXP nBurn, SEXP nFilter,SEXP nSubjects,SEXP nPredictSubjects){
    
    string fName = Rcpp::as<string>(fileName);

    int nS = Rcpp::as<int>(nSweeps);
    int nB = Rcpp::as<int>(nBurn);
    int nF = Rcpp::as<int>(nFilter);
    int nSj = Rcpp::as<int>(nSubjects);
    int nPSj = Rcpp::as<int>(nPredictSubjects);

    // Calculate how many samples we will have
    int nLines = 1 + (nS+nB)/nF;
    // Calculate when the burn in ends
    int firstLine = 2+nB/nF;

    vector<double> disSimMat((nSj*(nSj-1))/2+nPSj*nSj,1.0);
    vector<vector<double> > disSimMatPerSweep(1+nLines-firstLine);
    // Open the file with sample in
    ifstream zFile;
    zFile.open(fName.c_str());

    // The denominator for the contributions
    double denom = 1+(double)nLines-(double)firstLine;

    // A temporary vector for storing the cluster from each sweep
    vector<int> clusterData(nSj+nPSj);
    for(int k=1;k<=nLines;k++){
    	if(k<firstLine){
        	// Ignore the burn in
    		for(int i=0;i<nSj+nPSj;i++){
    			int tmp=0;
    			zFile >> tmp;
    		}
     	}else{
    		if((1+k-firstLine)==1||(1+k-firstLine)%100==0){
				std::cout << 1+k-firstLine << " samples out of " << 1+nLines-firstLine << std::endl;
			}
    		disSimMatPerSweep[k-firstLine+1].resize((nSj*(nSj-1))/2,1.0);
			for(int i=0;i<nSj+nPSj;i++){
    			// Fill up the cluster data for this sweep
    			zFile >> clusterData[i];
    		}


    		// Now we need to populate the dissimilarity matrix
    		int r=0;
			for(int i=0;i<nSj-1;i++){
    			for(int ii=i+1;ii<nSj;ii++){
    				if(clusterData[i]==clusterData[ii]){
        				disSimMat[r]-=1.0/denom;
        				disSimMatPerSweep[k-firstLine+1][r]-=1.0;
        			}
    				r++;
    			}
    		}

			for(int i=0;i<nPSj;i++){
				for(int ii=0;ii<nSj;ii++){
					if(clusterData[nSj+i]==clusterData[ii]){
						disSimMat[r]-=1.0/denom;
					}
					r++;
				}
			}
    	}
    }
    zFile.close();

    int minIndex=0;
    double currMinSum=nSj*(nSj-1)/2.0;
    for(int k=firstLine;k<=nLines;k++){
    	double tmpSum=0.0;
    	for(unsigned int r=0;r<nSj*(nSj-1)/2;r++){
    		tmpSum+=pow(disSimMat[r]-disSimMatPerSweep[k-firstLine+1][r],2.0);
    	}
    	if(tmpSum<currMinSum){
    		minIndex=k;
    		currMinSum=tmpSum;
    	}
    }


    return Rcpp::List::create(Rcpp::Named("lsOptSweep")=minIndex,
    						  Rcpp::Named("disSimMat")=Rcpp::wrap<vector<double> >(disSimMat));
}
