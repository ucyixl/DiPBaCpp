/// \file chain.h
/// \author David Hastie
/// \date 30 Sep 2010
/// \brief Header file defining classes for Markov chain Monte Carlo chains.

#ifndef CHAIN_H_
#define CHAIN_H_

// Standard includes
#include<fstream>
#include<string>
#include<vector>

// Custom includes
#include "MCMCpp/state.h"

using std::vector;
using std::cout;
using std::ostream;

/// \class mcmcChain chain.h "MCMC/chain.h"
/// \brief Class for Markov chain Monte Carlo chains
template<class modelParamType> class mcmcChain{

	public:
		/// \brief Default constructor
		mcmcChain() {};

		/// \brief Explicit constructor
		/// \param[in] currentState The initial state of the chain
		mcmcChain(const mcmcState<modelParamType>& currentState) :
			_currentState(currentState) {};

		/// \brief Destructor
		~mcmcChain(){};

		/// \brief Member function for setting the current state of the chain
		/// \param[in] x The state to be made the current state
		void currentState(const mcmcState<modelParamType>& x){
			_currentState = x;
		}

		/// \brief Member function to get the current state of the chain
		/// \return The current state of the chain
		const mcmcState<modelParamType>& currentState() const{
			return _currentState;
		}

		/// \brief Member function to get the current state of the chain
		/// \return The current state of the chain
		mcmcState<modelParamType>& currentState(){
			return _currentState;
		}


	private:
		/// \brief Current value of the chain
		mcmcState<modelParamType> _currentState;

};


#endif /* CHAIN_H_ */
