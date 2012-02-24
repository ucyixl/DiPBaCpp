/// \file DiPBaCOptions.h
/// \author David Hastie
/// \brief Header file for options class for DiPBaCpp

#ifndef DIPBACOPTIONS_H_
#define DIPBACOPTIONS_H_

// Standard includes
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<ctime>

using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;
using std::time;

/// \class diPBaCOptions DiPBaCOptions.h "DiPBaCOptions.h"
/// \brief A class for DiPBaC options
class diPBaCOptions{

	public:
		/// \brief Default constructor
		diPBaCOptions() {
			/* ---------- Default Values -------------------*/
			// Default values for the command line arguments
			// File names
			_inFileName = "input.txt";
			_outFileStem  = "output";
			_hyperParamFileName = "";
			_predictFileName = "";

			// MCMC sampler variables
			_nSweeps=10000;
			_nBurn=1000;
			_nProgress=500;
			_nFilter=1;
			// The default of 0 initial clusters, means a random number between
			// 5 and 15 is used at initialisation
			_nClusInit=0;

			// Random number seed
			_seed=(long) time(0);

			// Profile regression variables
			_outcomeType="Bernoulli";
			_covariateType="Discrete";
			_includeResponse = true;
			_responseExtraVar = false;
			_doPrediction = false;
			_varSelectType ="None";
			_fixedAlpha=-1;
		};

		/// \brief Default destructor
		~diPBaCOptions(){};

		/// \brief Return the number of sweeps
		unsigned int nSweeps() const{
			return _nSweeps;
		}

		/// \brief Set the number of sweeps
		void nSweeps(const unsigned int& nSw){
			_nSweeps=nSw;
		}

		/// \brief Return the number of burn in sweeps
		unsigned int nBurn() const{
			return _nBurn;
		}

		/// \brief Set the number of burn in sweeps
		void nBurn(const unsigned int& nB){
			_nBurn=nB;
		}

		/// \brief Return the number of filter sweeps
		unsigned int nFilter() const{
			return _nFilter;
		}

		/// \brief Set the number of filter sweeps
		void nFilter(const unsigned int& nF){
			_nFilter=nF;
		}

		/// \brief Return the number of progress iterations
		unsigned int nProgress() const{
			return _nProgress;
		}

		/// \brief Set the number of sweeps
		void nProgress(const unsigned int& nPr){
			_nProgress=nPr;
		}

		/// \brief Return the number of progress iterations
		unsigned int nClusInit() const{
			return _nClusInit;
		}

		/// \brief Set the number of sweeps
		void nClusInit(const unsigned int& nCl){
			_nClusInit=nCl;
		}

		/// \brief Return the random number seed
		long seed() const{
			return _seed;
		}

		/// \brief Set the random number seed
		void seed(const long& rndSeed){
			long rndSeedNew = rndSeed;
			if(rndSeedNew==0){
				rndSeedNew=(long) time(0);
			}
			_seed=rndSeedNew;
		}

		/// \brief Return the input file name
		string inFileName() const{
			return _inFileName;
		}

		/// \brief Set the input file name
		void inFileName(const string& inFile){
			_inFileName=inFile;
		}

		/// \brief Return the output file stem
		string outFileStem() const{
			return _outFileStem;
		}

		/// \brief Set the output file stem
		void outFileStem(const string& outStem){
			_outFileStem=outStem;
		}

		/// \brief Return the input file name
		string hyperParamFileName() const{
			return _hyperParamFileName;
		}

		/// \brief Set the input file name
		void hyperParamFileName(const string& hyperFile){
			_hyperParamFileName=hyperFile;
		}

		/// \brief Return the prediction file name
		string predictFileName() const{
			return _predictFileName;
		}

		/// \brief Set the prediction file name
		void predictFileName(const string& predictFile){
			_predictFileName=predictFile;
			_doPrediction=true;
		}

		/// \brief Return the outcome type
		string outcomeType() const{
			return _outcomeType;
		}

		/// \brief Set the outcome type
		void outcomeType(const string& outType){
			_outcomeType=outType;
		}

		/// \brief Return the covariate type
		string covariateType() const{
			return _covariateType;
		}

		/// \brief Set the covariate type
		void covariateType(const string& covType){
			_covariateType=covType;
		}

		/// \brief Return whether we are including response
		bool includeResponse() const{
			return _includeResponse;
		}

		/// \brief Set whether we are including the response
		void includeResponse(const bool& incResp){
			_includeResponse=incResp;
		}

		/// \brief Return the fixed value of alpha (negative value means not fixed)
		double fixedAlpha() const{
			return _fixedAlpha;
		}

		/// \brief Return the fixed value of alpha (negative value means not fixed)
		void fixedAlpha(const double& alphaVal){
			_fixedAlpha=alphaVal;
		}

		/// \brief Return whether we are including response
		bool doPrediction() const{
			return _doPrediction;
		}

		/// \brief Return whether there is response extra variation
		bool responseExtraVar() const{
			return _responseExtraVar;
		}

		/// \brief Set whether there is response extra variation
		void responseExtraVar(const bool& respExVar){
			_responseExtraVar=respExVar;
		}

		/// \brief Return the variable selection type
		string varSelectType() const{
			return _varSelectType;
		}

		/// \brief Set the variable selection type
		void varSelectType(const string& varSelType){
			_varSelectType=varSelType;
		}

		/// \brief Return whether we are resuming the run
		bool resumeRun() const{
			return _resumeRun;
		}

		/// \brief Set whether we are resuming the run
		void resumeRun(const bool& resume){
			_resumeRun=resume;
		}

		// Copy operator
		diPBaCOptions& operator=(const diPBaCOptions& options){

			_inFileName = options.inFileName();
			_outFileStem = options.outFileStem();
			_hyperParamFileName = options.hyperParamFileName();
			_predictFileName = options.predictFileName();
			_nSweeps = options.nSweeps();
			_nBurn = options.nBurn();
			_nFilter=options.nFilter();
			_nProgress=options.nProgress();
			_nClusInit=options.nClusInit();
			_seed=options.seed();
			_outcomeType=options.outcomeType();
			_covariateType=options.covariateType();
			_includeResponse=options.includeResponse();
			_fixedAlpha=options.fixedAlpha();
			_doPrediction=options.doPrediction();
			_responseExtraVar=options.responseExtraVar();
			_varSelectType=options.varSelectType();
			_resumeRun=options.resumeRun();
			return *this;
		}


	private:
		// The input file name
		string _inFileName;
		// The output file stem
		string _outFileStem;
		// The hyper parameter file
		string _hyperParamFileName;
		// The predictions file
		string _predictFileName;
		// The number of sweeps
		unsigned int _nSweeps;
		// The length of burn in
		unsigned int _nBurn;
		// How often the output should be filtered for reporting
		unsigned int _nFilter;
		// The number of iterations to print how sampler is progressing
		unsigned int _nProgress;
		// The number of initial clusters to initialise into
		unsigned int _nClusInit;
		// The random number seed
		long _seed;
		// The model for the outcome
		string _outcomeType;
		// The model for the covariates
		string _covariateType;
		// This notes whether we are including the response
		bool _includeResponse;
		// This has a fixed value of alpha (if negative we update alpha)
		double _fixedAlpha;
		// This notes whether we are also doing predictions
		bool _doPrediction;
		// This notes whether we have extra variation in the response
		bool _responseExtraVar;
		// The type of extra variation in the response
		string _varSelectType;
		// This notes whether we are resuming the run from a previous run
		bool _resumeRun;

};

#endif // DIPBACOPTIONS_H_
