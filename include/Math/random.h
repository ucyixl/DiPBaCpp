/// \file random.h
/// \author David Hastie
/// \date 1 Aug 2010
/// \brief Header file for generating random nos

#ifndef RANDOM_H_
#define RANDOM_H_

#include<string>
#include<iostream>
#include<limits>

#include<boost/random.hpp>
#include<boost/math/distributions/normal.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

using namespace Eigen;

using std::vector;
using std::string;

// Define the underlying random number generator
// mt19937 is the mersenne twister which is good for U(0,1) in up to 623 dimensions
typedef boost::random::mt19937 baseGeneratorType;

// Define the uniform random number generator
typedef boost::random::uniform_real_distribution<> randomUniform;
// Define the normal random number generator
typedef boost::random::normal_distribution<> randomNormal;
// Define the gamma random number generator
typedef boost::random::gamma_distribution<> randomGamma;
// Define the students t random number generator
typedef boost::random::students_t_distribution<> randomStudentsT;

// Univariate distributions

double betaRand(baseGeneratorType& rndGenerator,const double& a,const double& b){

	// Method taken from http://statprob.com/?op=getobj&from=objects&id=205
	randomGamma gammaRandA(a,1.0);
	randomGamma gammaRandB(b,1.0);

	double x1,x2;
	x1=gammaRandA(rndGenerator);
	x2=gammaRandB(rndGenerator);
	double out=x1/(x1+x2);
	return out;

}

double truncNormalRand(baseGeneratorType& rndGenerator,const double& mean,
					const double& stdDev,const string& distType,
					const double& lower,const double & upper){

	// Method as follows.
	// 1. Calculate the cdf values of lower and upper.
	// 2. Sample a prob uniformly between these values
	// 3. Use inverse cdf to transform p into sample
	// Note if distType="U" the distribution is upper truncated,
	// if distType="B" it is upper and lower truncated, and if it is "L"
	// it is lower truncated

	// Step 1.
	boost::math::normal_distribution<double> normDist(mean,stdDev);
	double pLower,pUpper,pSample;
	if(distType.compare("U")==0){
		pLower=0.0000000001;
		pUpper=cdf(normDist,upper);
	}else if(distType.compare("L")==0){
		pLower=cdf(normDist,lower);
		pUpper=1.0-0.0000000001;
	}else{
		pLower=cdf(normDist,lower);
		pUpper=cdf(normDist,upper);
	}

	// Step 2.
	// Define a uniform random number generator
	randomUniform unifRand(pLower,pUpper);
	pSample = unifRand(rndGenerator);
	// Step 3.
	return quantile(normDist,pSample);

}

// Multivariate distributions
vector<double> dirichletRand(baseGeneratorType& rndGenerator,const vector<double>& alpha){

	// We use the method from p. 482 of "Bayesian Data Analysis", Gelman et al.
	// Length of the working vector
	unsigned int n = alpha.size();

	vector<double> outVec(n);
	double sumVal = 0.0;
	for(unsigned int i=0;i<n;i++){
		randomGamma gammaRand(alpha[i],1.0);
		outVec[i]=gammaRand(rndGenerator);
		sumVal+=outVec[i];
	}

	for(unsigned int i=0;i<n;i++){
		outVec[i]/=sumVal;
	}

	return outVec;
}


VectorXd multivarNormalRand(baseGeneratorType& rndGenerator,const VectorXd& meanVec,const MatrixXd& covMat){

	unsigned int dimV = meanVec.size();

	// Create a normal random generator
	randomNormal normRand(0,1);

	VectorXd V(dimV);
	for(unsigned int i=0;i<dimV;i++){
		V(i)=normRand(rndGenerator);
	}

	// Cholesky decomposition
	LLT<MatrixXd> llt;
	MatrixXd B = (llt.compute(covMat)).matrixL();
	V = meanVec + B*V;
	return V;

}

MatrixXd invWishartRand(baseGeneratorType& rndGenerator,const MatrixXd& R,const int& m){

	// We generate a Wishart(invR,m) matrix variate then invert
	MatrixXd invR = R.inverse();

	MatrixXd invSample = wishartRand(rndGenerator,invR,m);

	return(invSample.inverse());

}


MatrixXd wishartRand(baseGeneratorType& rndGenerator,const MatrixXd& R,const int& m){

	// Cholesky decomposition
	LLT<MatrixXd> llt;
	MatrixXd D = (llt.compute(R)).matrixL();

	// Create a normal random generator
	randomNormal normRand(0,1);

	// Create the matrix A (s.t AA'~Wish(I,m))
	unsigned int dimR = R.rows();
	MatrixXd A=MatrixXd::Zero(dimR,dimR);

	// A Chi-squared with n d.o.f. is a Gamma(n/2,2)
	// using the shape and scale parametrisation used by Boost
	for(unsigned int i=0;i<dimR;i++){
		for(unsigned int j=0;j<i;j++){
			A(i,j)=normRand(rndGenerator);
		}
		randomGamma gammaRand((double)(m-i)/2.0,2.0);

		A(i,i)=sqrt(gammaRand(rndGenerator));
	}

	// Compute DA
	// Note DAA'D' = (DA)(DA)' ~ Wish(R,m)
	MatrixXd DA = D*A;

	// Return matrix square
	return DA*(DA.transpose());

}


#endif /* RANDOM_H_ */
