/// \file state.h
/// \author David Hastie
/// \date 26 Jul 2010
/// \brief Header file defining classes for Markov chain Monte Carlo state.

#ifndef STATE_H_
#define STATE_H_

// Standard includes
#include<vector>

// Use the standard library vectors
using std::vector;

/// \class mcmcState state.h "MCMC/state.h"
/// \brief Class for Markov chain Monte Carlo states
template<class modelParamType> class mcmcState {

	public:
		/// \brief Default constructor
		mcmcState() : _parameters(), _logPosterior(0.0), _logLikelihood(0.0),
			_logPrior(0.0) {};

		/// \brief Explicit constructor
		mcmcState(const modelParamType& parameters) : _logPosterior(0.0),
				_logLikelihood(0.0), _logPrior(0.0) {
			_parameters=parameters;
		};

		/// \brief Explicit constructor
		mcmcState(const modelParamType& parameters,double logPosterior,
				double logLikelihood,double logPrior) : _logPosterior(logPosterior),
				_logLikelihood(logLikelihood),_logPrior(logPrior){
			_parameters=parameters;
		};

		/// \brief Explicit constructor
		mcmcState(const modelParamType& parameters,const vector<double>& logPostVec) :
				_logPosterior(logPostVec[0]),_logLikelihood(logPostVec[1]),
				_logPrior(logPostVec[2]){
					_parameters=parameters;
		};

		const modelParamType& parameters() const{
			return _parameters;
		}

		modelParamType& parameters(){
			return _parameters;
		}


		// Destructor
		~mcmcState(){};


		/// \brief Member function to return the log posterior of the mcmc state
		/// \return The log posterior of the state
		double logPosterior() const {return _logPosterior;}

		/// \brief Member function to set the log posterior of the mcmc state
		void logPosterior(const double& logPost) {_logPosterior=logPost;}

		/// \brief Member function to return the logLikelihood of the mcmc state
		/// \return The log likelihood of the state
		double logLikelihood() const {return _logLikelihood;}

		/// \brief Member function to set the log likelihood of the mcmc state
		void logLikelihood(const double& logLik) {_logLikelihood=logLik;}

		/// \brief Member function to return the log prior of the mcmc state
		/// \return The log prior of the mcmc State
		double logPrior() const {return _logPrior;}

		/// \brief Member function to set the log Prior of the mcmc state
		void logPrior(const double& logPri) {_logPrior=logPri;}

		/// \brief Member function to set the log posterior, log likelihood
		/// and log prior of the mcmcState
		void logPosterior(const vector<double>& logPostVec){
			_logPosterior = logPostVec[0];
			_logLikelihood = logPostVec[1];
			_logPrior = logPostVec[2];
		}

		/// \brief Copy operator
		mcmcState& operator=(const mcmcState& x){
			_parameters=x._parameters;
			_logPosterior=x._logPosterior;
			_logLikelihood=x._logLikelihood;
			_logPrior=x._logPrior;
			return *this;
		}


	private:
		/// \brief The underlying MCMC state
		modelParamType _parameters;

		/// \brief The log posterior of the MCMC state
		double _logPosterior;

		/// \brief The log likelihood of the MCMC state
		double _logLikelihood;

		/// \brief The log prior of the MCMC state
		double _logPrior;

};



#endif /* STATE_H_ */
