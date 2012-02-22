/// \file proposal.h
/// \author David Hastie
/// \date 30 Sep 2010
/// \brief Header file defining classes for Markov chain Monte Carlo proposal.

#ifndef PROPOSAL_H_
#define PROPOSAL_H_

// Standard includes
#include<vector>
#include<string>

// Custom includes
#include<Math/random.h>
#include<MCMC/model.h>
#include<MCMC/chain.h>

using std::vector;
using std::string;

/// \class mcmcProposal proposal.h "MCMC/proposal.h"
/// \brief Class for Markov chain Monte Carlo proposal
/// \param[in] modelParamType The user defined type to contain model parameters
/// \param[in] dataType The user defined type to contain a data observation
///	(can be omitted for problems with no data)
template<class modelParamType,class optionType,class propParamType,class dataType = int> class mcmcProposal{

	public:
		/// \brief Default constructor
		mcmcProposal(){};

		/// \brief Destructor
		~mcmcProposal(){};

		/// \brief Member function to initialise the proposal
		/// \param[in] proposalParameters pointer to propParamType object
		/// containing initial parameter values
		/// \param[in] weight The probability of using this proposal on a
		/// sweep where it is attempted
		/// \param[in] frequency The proposal will be attempted after every
		/// "frequency" sweeps
		/// \note First the sampler checks whether sweep MOD frequency == 0 and if
		/// so the proposal is attempted with probability "weight". For example
		/// if frequency = 10 and weight=0.5, on average the proposal will be
		/// used every 20 sweeps.
		void initialiseProposal(const string& propName,
									propParamType* proposalParameters,
									const double& weight,
									const int& frequency,
									const int& firstSweep){
			_proposalName=propName;
			_proposalParameters = proposalParameters;
			_proposalWeight = weight;
			_proposalFrequency = frequency;
			_proposalFirstSweep = firstSweep;
			_nTry=0;
			_nAccept=0;
		}

		/// \brief Member function to set the function for calculating proposed
		/// new model parameters
		/// \param[in] *f Pointer to function to calculate proposed new model
		/// parameters
		void updateFunction(void (*f)(mcmcChain<modelParamType>&,
										unsigned int&,
										unsigned int&,
										const mcmcModel<modelParamType,optionType,dataType>&,
										propParamType&,
										baseGeneratorType&)){
			_updateParameters=f;
		}


		/// \brief Member function to return the probability of using this
		/// proposal given it is attempted
		/// \return The probability of using the proposal
		double proposalWeight(){
			return _proposalWeight;
		}

		/// \brief Member function to return the name of this proposal
		/// \return The name of the proposal
		string proposalName(){
			return _proposalName;
		}

		/// \brief Member function to return the frequency (in sweeps)
		/// that this proposal should be attempted
		/// \return The frequency in sweeps that the proposal should be attempted
		unsigned int proposalFrequency(){
			return _proposalFrequency;
		}

		/// \brief Member function to return the first sweep
		/// that this proposal should be attempted
		/// \return The first sweep that the proposal should be attempted
		unsigned int proposalFirstSweep(){
			return _proposalFirstSweep;
		}

		/// \brief Member function to propose new model functions
		/// \param[in,out] chain A reference to the mcmcChain
		/// \param[in] rndGenerator A boost random number generator object
		void updateParameters(mcmcChain<modelParamType>& chain,
							  const mcmcModel<modelParamType,optionType,dataType>& model,
										baseGeneratorType& rndGenerator){

			// Use the user function to propose some new parameters
			(*_updateParameters)(chain,
								_nTry,
								_nAccept,
								model,
								*_proposalParameters,
								rndGenerator);

		}

		/// \brief Member function to return the acceptance rate
		double acceptanceRate() const{
			if(_nTry>0){
				return (double)_nAccept / (double)_nTry;
			}else{
				return 0.0;
			}
		}

		/// \brief Member function to return the number of tries of the proposal
		unsigned int nTry() const{
			return _nTry;
		}

		/// \brief Member function to return the number of acceptances of proposal
		unsigned int nAccept() const{
			return _nAccept;
		}

	private:
		string _proposalName;

		/// \brief The parameter used for the proposal
		propParamType* _proposalParameters;

		/// \var _nTry
		/// \brief The number of times the move has been tried
		/// \var _nAccept
		/// \brief The number of times the move has been accepted
		unsigned int _nTry;
		unsigned int _nAccept;

		/// \brief Pointer to user function to propose new model parameters
		void (*_updateParameters)(mcmcChain<modelParamType>&,
								unsigned int &,
								unsigned int&,
								const mcmcModel<modelParamType,optionType,dataType>&,
								propParamType&,
								baseGeneratorType&);

		/// \brief The propbability that the proposal is used given that it is attempted
		double _proposalWeight;

		/// \brief The frequency with which the proposal is attempted
		unsigned int _proposalFrequency;

		/// \brief The first sweep at which the proposal should be attempted
		unsigned int _proposalFirstSweep;

};

#endif /* PROPOSAL_H_ */
