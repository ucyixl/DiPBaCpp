/// \file DiPBaC.cpp
/// \author David Hastie
/// \brief Main file for running DiPBaC

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
#include "DiPBaCOptions.h"
#include "DiPBaCModel.h"
#include "DiPBaCData.h"
#include "DiPBaCProposals.h"
#include "DiPBaCIO.h"

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
	diPBaCOptions options = processCommandLine(argc,argv);
	
	/* ---------- Set up the sampler object--------*/
	// Initialise the sampler object
	mcmcSampler<diPBaCParams,diPBaCOptions,
				diPBaCPropParams,diPBaCData> diPBaCSampler;

	// Set the options
	diPBaCSampler.options(options);

	// Set the model
	diPBaCSampler.model(&importDiPBaCData,&initialiseDiPBaC,
							&diPBaCLogPost,true);

	// Set the missing data function
	diPBaCSampler.updateMissingDataFn(&updateMissingDiPBaCData);

	// Add the function for writing output
	diPBaCSampler.userOutputFn(&writeDiPBaCOutput);

	// Seed the random number generator
	diPBaCSampler.seedGenerator(options.seed());

	// Set the sampler specific variables
	diPBaCSampler.nSweeps(options.nSweeps());
	diPBaCSampler.nBurn(options.nBurn());
	diPBaCSampler.nFilter(options.nFilter());
	diPBaCSampler.nProgress(options.nProgress());
	diPBaCSampler.reportBurnIn(true);

	/* ---------- Read in the data -------- */
	diPBaCSampler.model().dataset().outcomeType(options.outcomeType());
	diPBaCSampler.model().dataset().covariateType(options.covariateType());
	diPBaCSampler.importData(options.inFileName(),options.predictFileName());
	diPBaCData dataset = diPBaCSampler.model().dataset();

	/* ---------- Add the proposals -------- */

	// Set the proposal parameters
	diPBaCPropParams proposalParams(options.nSweeps(),dataset.nCovariates(),
										dataset.nConfounders());
	diPBaCSampler.proposalParams(proposalParams);

	// The gibbs update for the active V
	diPBaCSampler.addProposal("gibbsForVActive",1.0,1,1,&gibbsForVActive);


	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		if(dataset.anyCategorical()){
			//  Update for the active phi parameters
			diPBaCSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);
		}

		if(dataset.anyOrdinal()){
			// Adaptive MH for active delta
			diPBaCSampler.addProposal("metropolisHastingsForDeltaActive",1.0,1,1,&metropolisHastingsForDeltaActive);
		}

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		diPBaCSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);

		// Update for the active Sigma parameters
		diPBaCSampler.addProposal("gibbsForTauActive",1.0,1,1,&gibbsForTauActive);

	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			diPBaCSampler.addProposal("gibbsForGammaActive",1.0,1,firstSweep,&gibbsForGammaActive);
		}

	}

	if(options.includeResponse()){
		// The Metropolis Hastings update for the active theta
		diPBaCSampler.addProposal("metropolisHastingsForThetaActive",1.0,1,1,&metropolisHastingsForThetaActive);
	}

	// The Metropolis Hastings update for labels
	diPBaCSampler.addProposal("metropolisHastingsForLabels",1.0,1,1,&metropolisHastingsForLabels);

	// Gibbs for U
	diPBaCSampler.addProposal("gibbsForU",1.0,1,1,&gibbsForU);

	// The Metropolis Hastings update for alpha
	if(options.fixedAlpha()<0){
		diPBaCSampler.addProposal("metropolisHastingsForAlpha",1.0,1,1,&metropolisHastingsForAlpha);
	}

	// The gibbs update for the inactive V
	diPBaCSampler.addProposal("gibbsForVInActive",1.0,1,1,&gibbsForVInActive);

	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		if(dataset.anyCategorical()){
			//  Update for the inactive phi parameters
			diPBaCSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);
		}

		if(dataset.anyOrdinal()){
			// Adaptive MH for inactive delta
			diPBaCSampler.addProposal("metropolisHastingsForDeltaInActive",1.0,1,1,&metropolisHastingsForDeltaInActive);
		}

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		diPBaCSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);

		// Update for the active Sigma parameters
		diPBaCSampler.addProposal("gibbsForTauInActive",1.0,1,1,&gibbsForTauInActive);
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			diPBaCSampler.addProposal("gibbsForGammaInActive",1.0,1,firstSweep,&gibbsForGammaInActive);
		}

	}

	if(options.includeResponse()){
		// The Metropolis Hastings update for the inactive theta
		diPBaCSampler.addProposal("gibbsForThetaInActive",1.0,1,1,&gibbsForThetaInActive);
	}

	if(options.includeResponse()){
		// Adaptive MH for beta
		if(dataset.nConfounders()>0){
			diPBaCSampler.addProposal("metropolisHastingsForBeta",1.0,1,1,&metropolisHastingsForBeta);
		}

		if(options.responseExtraVar()){
			// Adaptive MH for lambda
			diPBaCSampler.addProposal("metropolisHastingsForLambda",1.0,1,1,&metropolisHastingsForLambda);

			// Gibbs for tauEpsilon
			diPBaCSampler.addProposal("gibbsForTauEpsilon",1.0,1,1,&gibbsForTauEpsilon);
		}
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		// Metropolis Hastings for joint update of rho and omega
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);

		diPBaCSampler.addProposal("metropolisHastingsForRhoOmega",1.0,1,firstSweep,&metropolisHastingsForRhoOmega);
	}


	// Gibbs update for the allocation parameters
	diPBaCSampler.addProposal("gibbsForZ",1.0,1,1,&gibbsForZ);


	/* ---------- Initialise the output files -----*/
	diPBaCSampler.initialiseOutputFiles(options.outFileStem());

	/* ---------- Write the log file ------------- */
	// The standard log file
	diPBaCSampler.writeLogFile();

	/* ---------- Initialise the chain ---- */
	diPBaCSampler.initialiseChain();
	diPBaCHyperParams hyperParams = diPBaCSampler.chain().currentState().parameters().hyperParams();
	/* ---------- Run the sampler --------- */
	// Note: in this function the output gets written
	diPBaCSampler.run();

	/* -- End the clock time and write the full run details to log file --*/
	currTime = time(NULL);
    double timeInSecs=(double)currTime-(double)beginTime;
	string tmpStr = storeLogFileData(options,dataset,hyperParams,timeInSecs);
	diPBaCSampler.appendToLogFile(tmpStr);


	/* ---------- Clean Up ---------------- */
	diPBaCSampler.closeOutputFiles();
	return(0);

}
