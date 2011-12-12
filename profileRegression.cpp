/// \file profileRegression.cpp
/// \author David Hastie
/// \date 26 Jan 2011
/// \brief Main file for running profile regression

// Standard includes
#include <cmath>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <sstream>

// Custom includes
#include "MCMC/sampler.h"
#include "MCMC/model.h"
#include "MCMC/proposal.h"
#include "profileRegression/profileRegressionOptions.h"
#include "profileRegression/profileRegressionModel.h"
#include "profileRegression/profileRegressionData.h"
#include "profileRegression/profileRegressionProposals.h"
#include "profileRegression/profileRegressionIO.h"

using std::vector;
using std::cout;
using std::endl;
using std::ostringstream;
using std::time;

int main(int argc, char*  argv[]){

	/* ---------- Start the timer ------------------*/
    time_t beginTime,currTime;
	beginTime = time(NULL);

	/* -----------Process the command line ---------*/
	profRegrOptions options = processCommandLine(argc,argv);
	
	/* ---------- Set up the sampler object--------*/
	// Initialise the sampler object
	mcmcSampler<profRegrParams,profRegrOptions,
				profRegrPropParams,profRegrData> profRegrSampler;

	// Set the options
	profRegrSampler.options(options);

	// Set the model
	profRegrSampler.model(&importProfRegrData,&initialiseProfRegr,
							&profRegrLogPost,true);

	// Set the missing data function
	profRegrSampler.updateMissingDataFn(&updateMissingProfRegrData);

	// Add the function for writing output
	profRegrSampler.userOutputFn(&writeProfRegrOutput);

	// Seed the random number generator
	profRegrSampler.seedGenerator(options.seed());

	// Set the sampler specific variables
	profRegrSampler.nSweeps(options.nSweeps());
	profRegrSampler.nBurn(options.nBurn());
	profRegrSampler.nFilter(options.nFilter());
	profRegrSampler.nProgress(options.nProgress());
	profRegrSampler.reportBurnIn(true);

	/* ---------- Read in the data -------- */
	profRegrSampler.model().dataset().outcomeType(options.outcomeType());
	profRegrSampler.model().dataset().covariateType(options.covariateType());
	profRegrSampler.importData(options.inFileName(),options.predictFileName());
	profRegrData dataset = profRegrSampler.model().dataset();

	/* ---------- Add the proposals -------- */

	// Set the proposal parameters
	profRegrPropParams proposalParams(options.nSweeps(),dataset.nCovariates(),
										dataset.nConfounders());
	profRegrSampler.proposalParams(proposalParams);

	// The gibbs update for the active V
	profRegrSampler.addProposal("gibbsForVActive",1.0,1,1,&gibbsForVActive);


	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		if(dataset.anyCategorical()){
			//  Update for the active phi parameters
			profRegrSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);
		}

		if(dataset.anyOrdinal()){
			// Adaptive MH for active delta
			profRegrSampler.addProposal("metropolisHastingsForDeltaActive",1.0,1,1,&metropolisHastingsForDeltaActive);
		}

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		profRegrSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);

		// Update for the active Sigma parameters
		profRegrSampler.addProposal("gibbsForTauActive",1.0,1,1,&gibbsForTauActive);

	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			profRegrSampler.addProposal("gibbsForGammaActive",1.0,1,firstSweep,&gibbsForGammaActive);
		}

	}

	if(options.includeResponse()){
		// The Metropolis Hastings update for the active theta
		profRegrSampler.addProposal("metropolisHastingsForThetaActive",1.0,1,1,&metropolisHastingsForThetaActive);
	}

	// The Metropolis Hastings update for labels
	profRegrSampler.addProposal("metropolisHastingsForLabels",1.0,1,1,&metropolisHastingsForLabels);

	// Gibbs for U
	profRegrSampler.addProposal("gibbsForU",1.0,1,1,&gibbsForU);

	// The Metropolis Hastings update for alpha
	if(options.fixedAlpha()<0){
		profRegrSampler.addProposal("metropolisHastingsForAlpha",1.0,1,1,&metropolisHastingsForAlpha);
	}

	// The gibbs update for the inactive V
	profRegrSampler.addProposal("gibbsForVInActive",1.0,1,1,&gibbsForVInActive);

	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		if(dataset.anyCategorical()){
			//  Update for the inactive phi parameters
			profRegrSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);
		}

		if(dataset.anyOrdinal()){
			// Adaptive MH for inactive delta
			profRegrSampler.addProposal("metropolisHastingsForDeltaInActive",1.0,1,1,&metropolisHastingsForDeltaInActive);
		}

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		profRegrSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);

		// Update for the active Sigma parameters
		profRegrSampler.addProposal("gibbsForTauInActive",1.0,1,1,&gibbsForTauInActive);
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			profRegrSampler.addProposal("gibbsForGammaInActive",1.0,1,firstSweep,&gibbsForGammaInActive);
		}

	}

	if(options.includeResponse()){
		// The Metropolis Hastings update for the inactive theta
		profRegrSampler.addProposal("gibbsForThetaInActive",1.0,1,1,&gibbsForThetaInActive);
	}

	if(options.includeResponse()){
		// Adaptive MH for beta
		if(dataset.nConfounders()>0){
			profRegrSampler.addProposal("metropolisHastingsForBeta",1.0,1,1,&metropolisHastingsForBeta);
		}

		if(options.responseExtraVar()){
			// Adaptive MH for lambda
			profRegrSampler.addProposal("metropolisHastingsForLambda",1.0,1,1,&metropolisHastingsForLambda);

			// Gibbs for tauEpsilon
			profRegrSampler.addProposal("gibbsForTauEpsilon",1.0,1,1,&gibbsForTauEpsilon);
		}
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		// Metropolis Hastings for joint update of rho and omega
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);

		profRegrSampler.addProposal("metropolisHastingsForRhoOmega",1.0,1,firstSweep,&metropolisHastingsForRhoOmega);
	}


	// Gibbs update for the allocation parameters
	profRegrSampler.addProposal("gibbsForZ",1.0,1,1,&gibbsForZ);


	/* ---------- Initialise the output files -----*/
	profRegrSampler.initialiseOutputFiles(options.outFileStem());

	/* ---------- Write the log file ------------- */
	// The standard log file
	profRegrSampler.writeLogFile();

	/* ---------- Initialise the chain ---- */
	profRegrSampler.initialiseChain();
	profRegrHyperParams hyperParams = profRegrSampler.chain().currentState().parameters().hyperParams();
	/* ---------- Run the sampler --------- */
	// Note: in this function the output gets written
	profRegrSampler.run();

	/* -- End the clock time and write the full run details to log file --*/
	currTime = time(NULL);
    double timeInSecs=(double)currTime-(double)beginTime;
	string tmpStr = storeLogFileData(options,dataset,hyperParams,timeInSecs);
	profRegrSampler.appendToLogFile(tmpStr);


	/* ---------- Clean Up ---------------- */
	profRegrSampler.closeOutputFiles();
	return(0);

}

