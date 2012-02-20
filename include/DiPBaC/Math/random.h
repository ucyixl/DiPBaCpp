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
#include <boost/math/constants/constants.hpp>

#include<Eigen/Dense>

#include "MCMCpp/random.h"

using namespace Eigen;
using namespace boost::math::constants;

using std::vector;
using std::string;
using std::cout;
using std::endl;

using boost::math::normal_distribution;


// Define the uniform random number generator
typedef boost::variate_generator<baseGeneratorType&,boost::uniform_real<double> > randomUniform;
// Define the normal random number generator
typedef boost::variate_generator<baseGeneratorType&,boost::normal_distribution<double> > randomNormal;
// Define the gamma random number generator
typedef boost::variate_generator<baseGeneratorType&,boost::gamma_distribution<double> > randomGamma;

unsigned int sampleSingleCumulative(baseGeneratorType& rndGenerator,const vector<double>& cumProbVec){

	boost::uniform_real<double> unifDist(0,1);
	randomUniform unifRand(rndGenerator,unifDist);
	if(fabs(cumProbVec[cumProbVec.size()-1]-1)>0.000000000001){
		cout << "Cumulative sampling distribution does not sum to 1" << endl;
		exit(-1);
	}
	unsigned int out=0;
	double u = unifRand();
	for(unsigned int i=0;i<cumProbVec.size();i++){
		if(u<cumProbVec[i]){
			out=i;
			break;
		}
	}
	return out;
}

unsigned int sampleSingleNormalised(baseGeneratorType& rndGenerator,const vector<double>& probVec){

	unsigned int vecSize = probVec.size();
	vector<double> cumProbVec(vecSize);
	unsigned int out=0;
	if(vecSize>0){
		cumProbVec[0]=probVec[0];
		for(unsigned int i=1;i<vecSize;i++){
			cumProbVec[i]=cumProbVec[i-1]+probVec[i];
		}
		out = sampleSingleCumulative(rndGenerator,cumProbVec);
	}
	return out;
}


unsigned int sampleSingleUnnormalised(baseGeneratorType& rndGenerator,const vector<double>& propProbVec){

	unsigned int vecSize = propProbVec.size();
	vector<double> cumProbVec(vecSize);
	unsigned int out=0;
	if(vecSize>0){
		double sum = 0.0;
		for(unsigned int i=0;i<vecSize;i++){
			sum+=propProbVec[i];
		}

		cumProbVec[0]=propProbVec[0]/sum;
		for(unsigned int i=1;i<vecSize;i++){
			cumProbVec[i]=cumProbVec[i-1]+propProbVec[i]/sum;
		}
		out = sampleSingleCumulative(rndGenerator,cumProbVec);
	}
	return out;

}

unsigned int sampleSingleCumulative(const double& u,const vector<double>& cumProbVec){

	if(fabs(cumProbVec[cumProbVec.size()-1]-1)>0.000000000001){
		cout << "Cumulative sampling distribution does not sum to 1" << endl;
		exit(-1);
	}
	unsigned int out=0;
	for(unsigned int i=0;i<cumProbVec.size();i++){
		if(u<cumProbVec[i]){
			out=i;
			break;
		}
	}
	return out;
}

unsigned int sampleSingleNormalised(const double& u,const vector<double>& probVec){

	unsigned int vecSize = probVec.size();
	vector<double> cumProbVec(vecSize);
	unsigned int out=0;
	if(vecSize>0){
		cumProbVec[0]=probVec[0];
		for(unsigned int i=1;i<vecSize;i++){
			cumProbVec[i]=cumProbVec[i-1]+probVec[i];
		}
		out = sampleSingleCumulative(u,cumProbVec);
	}
	return out;
}


unsigned int sampleSingleUnnormalised(const double& u,const vector<double>& propProbVec){

	unsigned int vecSize = propProbVec.size();
	vector<double> cumProbVec(vecSize);
	unsigned int out=0;
	if(vecSize>0){
		double sum = 0.0;
		for(unsigned int i=0;i<vecSize;i++){
			sum+=propProbVec[i];
		}

		cumProbVec[0]=propProbVec[0]/sum;
		for(unsigned int i=1;i<vecSize;i++){
			cumProbVec[i]=cumProbVec[i-1]+propProbVec[i]/sum;
		}
		out = sampleSingleCumulative(u,cumProbVec);
	}
	return out;

}


double betaRand(baseGeneratorType& rndGenerator,const double& a,const double& b){

	// Method taken from http://statprob.com/?op=getobj&from=objects&id=205
	boost::gamma_distribution<double> gammaDistA(a);
	boost::gamma_distribution<double> gammaDistB(b);
	randomGamma gammaRandA(rndGenerator,gammaDistA);
	randomGamma gammaRandB(rndGenerator,gammaDistB);

	double x1,x2;
	x1=gammaRandA();
	x2=gammaRandB();
	double out=x1/(x1+x2);
	if(std::fabs(1.0 - out) < 5*std::numeric_limits<double>::epsilon()){
		out=1.0-5*std::numeric_limits<double>::epsilon();
	}
	return out;

}

double expRand(baseGeneratorType& rndGenerator,const double& beta){

	boost::uniform_real<double> unifDist(0,1);
	randomUniform unifRand(rndGenerator,unifDist);

	return -log(unifRand())/beta;

}

double studentsTRand(baseGeneratorType& rndGenerator,const unsigned int& dof){

	// Method from http://statprob.com/?op=getobj&from=objects&id=205
	// Generate a standard normal
	boost::normal_distribution<double> normDist(0,1);
	randomNormal normRand(rndGenerator,normDist);
	double x1=normRand();

	// Generate a Chi squared
	// Note: if X~ Gamma(a,1), then X/b~Gamma(a,b), using
	// the shape and inverse scale parameterisation
	// A Chi-squared with n d.o.f. is a Gamma(n/2,1/2)
	boost::gamma_distribution<double> gammaDist((double)dof/2.0);
	randomGamma gammaRand(rndGenerator,gammaDist);
	double x2 = 2*gammaRand();

	return x1/sqrt(x2/(double)dof);


}

vector<double> dirichletRand(baseGeneratorType& rndGenerator,const vector<double>& alpha){

	// We use the method from p. 482 of "Bayesian Data Analysis", Gelman et al.
	// Length of the working vector
	unsigned int n = alpha.size();

	vector<double> outVec(n);
	double sumVal = 0.0;
	for(unsigned int i=0;i<n;i++){
		boost::gamma_distribution<double> gammaDist(alpha[i]);
		randomGamma gammaRand(rndGenerator,gammaDist);
		outVec[i]=gammaRand();
		sumVal+=outVec[i];
	}

	for(unsigned int i=0;i<n;i++){
		outVec[i]/=sumVal;
	}

	return outVec;
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
	normal_distribution<double> normDist(mean,stdDev);
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
	boost::uniform_real<double> unifDist(pLower,pUpper);
	randomUniform unifRand(rndGenerator,unifDist);
	pSample = unifRand();
	// Step 3.
	return quantile(normDist,pSample);

}

VectorXd multivarNormalRand(baseGeneratorType& rndGenerator,const VectorXd& meanVec,const MatrixXd& covMat){

	unsigned int dimV = meanVec.size();

	// Create a normal random generator
	boost::normal_distribution<double> normDist(0,1);
	randomNormal normRand(rndGenerator,normDist);

	VectorXd V(dimV);
	for(unsigned int i=0;i<dimV;i++){
		V(i)=normRand();
	}

	// Cholesky decomposition
	LLT<MatrixXd> llt;
	MatrixXd B = (llt.compute(covMat)).matrixL();
	V = meanVec + B*V;
	return V;

}


MatrixXd wishartRand(baseGeneratorType& rndGenerator,const MatrixXd& R,const int& m){

	// Cholesky decomposition
	LLT<MatrixXd> llt;
	MatrixXd D = (llt.compute(R)).matrixL();

	// Create a normal random generator
	boost::normal_distribution<double> normDist(0,1);
	randomNormal normRand(rndGenerator,normDist);

	// Create the matrix A (s.t AA'~Wish(I,m))
	unsigned int dimR = R.rows();
	MatrixXd A=MatrixXd::Zero(dimR,dimR);

	// Note: if X~ Gamma(a,1), then X/b~Gamma(a,b), using
	// the shape and inverse scale parameterisation
	// A Chi-squared with n d.o.f. is a Gamma(n/2,1/2)
	for(unsigned int i=0;i<dimR;i++){
		for(unsigned int j=0;j<i;j++){
			A(i,j)=normRand();
		}
		boost::gamma_distribution<double> gammaDist((double)(m-i)/2.0);
		randomGamma gammaRand(rndGenerator,gammaDist);

		A(i,i)=sqrt(2*gammaRand());
	}

	// Compute DA
	// Note DAA'D' = (DA)(DA)' ~ Wish(R,m)
	MatrixXd DA = D*A;

	// Return matrix square
	return DA*(DA.transpose());

}

MatrixXd invWishartRand(baseGeneratorType& rndGenerator,const MatrixXd& R,const int& m){

	// We generate a Wishart(invR,m) matrix variate then invert
	MatrixXd invR = R.inverse();

	MatrixXd invSample = wishartRand(rndGenerator,invR,m);

	return(invSample.inverse());

}

#endif /* RANDOM_H_ */
