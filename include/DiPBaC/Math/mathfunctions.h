/// \file mathfunctions.h
/// \author David Hastie
/// \date 18 Oct 2010
/// \brief Header file to define distributions

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include<cmath>

#include<boost/math/special_functions/gamma.hpp>

using namespace boost::math::constants;

using boost::math::lgamma;

double logMultivarGammaFn(const double& x,const unsigned int& p){

	double out;
	out = 0.25*(double)(p*(p-1))*log(pi<double>());
	for(unsigned int i=1;i<=p;i++){
		out += lgamma(x+(1.0-(double)i)/2.0);
	}
	return out;
}

double logit(const double& lambda){
	return 1.0/(1.0+exp(-lambda));
}

#endif /*MATHFUNCTIONS_H_*/
