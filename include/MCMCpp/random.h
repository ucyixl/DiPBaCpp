/// \file random.h
/// \author David Hastie
/// \date 1 Aug 2010
/// \brief Header file for generating random nos

#ifndef MCMCPP_RANDOM_H_
#define MCMCPP_RANDOM_H_

#include<boost/random.hpp>

// Define the underlying random number generator
// mt19937 is the mersenne twister which is good for U(0,1) in up to 623 dimensions
typedef boost::mt19937 baseGeneratorType;

#endif /* MCMCPP_RANDOM_H_ */
