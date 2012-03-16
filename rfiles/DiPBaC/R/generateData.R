
# Generate simulated data for testing C++ DiPBaC
generateSampleDataFile<-function(clusterSummary,fileName){

	
	subjectsPerCluster<-clusterSummary$clusterSizes
	nSubjects<-sum(subjectsPerCluster)
   nCovariates<-clusterSummary$nCovariates
   covariateType<-clusterSummary$covariateType
   missingDataProb<-clusterSummary$missingDataProb
   nFixedEffects<-clusterSummary$nFixedEffects
   
	# Clustering covariates X
	X<-matrix(-999,nSubjects,nCovariates)
	
   k<-1
   
   # Loop over subjects
   for(i in 1:nSubjects){
      if(i<=subjectsPerCluster[k]){
         clusterData<-clusterSummary$clusterData[[k]]
      }else{
         clusterData<-clusterSummary$clusterData[[k+1]]         
         k<-k+1
         subjectsPerCluster[k]<-subjectsPerCluster[k]+subjectsPerCluster[k-1]
      }
      
      # Loop over covariates to generate the X data
      if(covariateType=='Discrete'){
         for(j in 1:nCovariates){
            if(i>1&&runif(1)<missingDataProb){
               X[i,j]<--999
            }else{
               u<-runif(1)
              	nCategories<-length(clusterData$covariateProbs[[j]])
               for(kk in 1:nCategories){
               	if(u<cumsum(clusterData$covariateProbs[[j]])[kk]){
                  	X[i,j]<-kk-1
                     break
                  }
               }
            }
         }
      }else if(covariateType=='Normal'){
         X[i,]<-clusterData$covariateMeans+t(chol(clusterData$covariateCovariance))%*%rnorm(nCovariates,0,1)   
         for(j in 1:nCovariates){
            if(i>1&&runif(1)<missingDataProb){
               X[i,j]<--999
            }
         }
      }
      
   }
	
	# Fixed Effects W
	if(nFixedEffects>0){
   	W<-matrix(rnorm(nSubjects*nFixedEffects,0,1),nSubjects,nFixedEffects)
   }else{
      W<-NULL
   }
		
	# Response Vector Y
   Y<-rep(0,nSubjects)
   outcomeType<-clusterSummary$outcomeType
   if(nFixedEffects>0){
      beta<-clusterSummary$fixedEffectsCoeffs   
   }
   
	if(outcomeType=='Poisson'){
      offset<-runif(nSubjects,clusterSummary$offsetLims[1],clusterSummary$offsetLims[2])
   }else{
      offset<-NULL
   }
   
	if(outcomeType=='Binomial'){
		nTrials<-sample(clusterSummary$nTrialsLims[1]:clusterSummary$nTrials[2],nSubjects,replace=T)
	}else{
		nTrials<-NULL
	}
	
	if(outcomeType=='Normal'){
		sigmaSqY=clusterSummary$sigmaSqY
	}
	
   subjectsPerCluster<-clusterSummary$clusterSizes
	k<-1
   # Loop over subjects
   for(i in 1:nSubjects){
            
      if(i<=subjectsPerCluster[k]){
         theta<-clusterSummary$clusterData[[k]]$theta
      }else{
         theta<-clusterSummary$clusterData[[k+1]]$theta         
         k<-k+1
         subjectsPerCluster[k]<-subjectsPerCluster[k]+subjectsPerCluster[k-1]
      }
		
		mu<-theta
      if(nFixedEffects>0){
         mu<-mu+sum(beta*W[i,])         
      }
      if(outcomeType=='Poisson'){
         mu<-mu+log(offset[i])
         Y[i]<-rpois(1,exp(mu))
      }else if(outcomeType=='Bernoulli'){
         p<-1/(1+exp(-mu))
         if(runif(1)<p){
            Y[i]<-1
         }else{
            Y[i]<-0
         }
      }else if(outcomeType=='Binomial'){
         p<-1/(1+exp(-mu))
			Y[i]<-sum(runif(nTrials[i])<p)
		}else if(outcomeType=='Normal'){
			Y[i]<-rnorm(1,mu,sqrt(sigmaSqY))
		}
      
   }

   # Write the output
   write(nSubjects,fileName,append=F)
   nCovariates<-clusterSummary$nCovariates
   nFixedEffects<-clusterSummary$nFixedEffects
   write(nCovariates,fileName,append=T)
   covNames<-paste('Variable',seq(1,nCovariates,1),sep="")
   write(covNames,fileName,append=T,ncolumns=1)
   write(nFixedEffects,fileName,append=T)
   if(nFixedEffects>0){
      fixEffNames<-paste('FixedEffects',seq(1,nFixedEffects,1),sep="")
      write(fixEffNames,fileName,append=T,ncolumns=1)
   }
   if(clusterSummary$covariateType=='Discrete'){
      write(clusterSummary$nCategories,fileName,ncolumns=nCovariates,append=T)
   }
   outData<-paste(Y,apply(X,1,paste,collapse=" "))
   if(nFixedEffects>0){
      outData<-paste(outData,apply(W,1,paste,collapse=" "))
   }
   if(clusterSummary$outcomeType=="Poisson"){
      outData<-paste(outData,offset)
   }
   if(clusterSummary$outcomeType=="Binomial"){
      outData<-paste(outData,nTrials)
   }
   write(outData,fileName,append=T,ncolumns=1)
}


############################
# Sample datasets


clusSummaryBernoulliDiscrete<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=5,
      'nCategories'=c(3,3,3,3,3),
      'nFixedEffects'=0,
      'missingDataProb'=0.001,
      'nClusters'=5,
		'clusterSizes'=c(100,200,300,200,100),
      'clusterData'=list(list('theta'=log(9),
                  'covariateProbs'=list(c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1))),
            list('theta'=log(2),
                  'covariateProbs'=list(c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1))),
            list('theta'=0,
                  'covariateProbs'=list(c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1))),
				list('theta'=log(1/2),
                  'covariateProbs'=list(c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.8,0.1),
                        c(0.1,0.1,0.8),
                        c(0.1,0.1,0.8))),
            
            list('theta'=log(1/9),
                  'covariateProbs'=list(c(0.1,0.1,0.8),
                        c(0.1,0.1,0.8),
                        c(0.1,0.1,0.8),
                        c(0.1,0.1,0.8),
                        c(0.1,0.1,0.8)))))

clusSummaryBernoulliDiscreteIdentical<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=5,
      'nCategories'=c(3,3,3,3,3),
      'nFixedEffects'=0,
      'missingDataProb'=0,
      'nClusters'=5,
      'clusterSizes'=c(100,200,300,200,100),
      'clusterData'=list(list('theta'=log(9),
                   'covariateProbs'=list(c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1))),
                    list('theta'=log(2),
                    'covariateProbs'=list(c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1))),
                    list('theta'=0,
                    'covariateProbs'=list(c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1))),
                    list('theta'=log(1/2),
                    'covariateProbs'=list(c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1))),                        
                    list('theta'=log(1/9),
                    'covariateProbs'=list(c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1),
                          c(0.8,0.1,0.1)))))

clusSummaryBernoulliDiscreteIdentical2<-list(
                    'outcomeType'='Bernoulli',
                    'covariateType'='Discrete',
                    'nCovariates'=2,
                    'nCategories'=c(3,3),
                    'nFixedEffects'=0,
                    'missingDataProb'=0,
                    'nClusters'=2,
                    'clusterSizes'=c(400,400),
                    'clusterData'=list(list('theta'=log(9),
                                'covariateProbs'=list(c(0.8,0.1,0.1),
                                      c(0.8,0.1,0.1))),
                          list('theta'=log(1/9),
                                'covariateProbs'=list(c(0.8,0.1,0.1),
                                      c(0.8,0.1,0.1)))))        
              
clusSummaryPoissonDiscrete<-list(
      'outcomeType'='Poisson',
      'covariateType'='Discrete',
      'nCovariates'=5,
      'nCategories'=c(2,2,3,3,4),
      'nFixedEffects'=1,
      'fixedEffectsCoeffs'=c(0.01),
      'offsetLims'=c(0.9,1.1),
      'missingDataProb'=0.001,
      'nClusters'=5,
		'clusterSizes'=c(150,250,250,250,150),
      'clusterData'=list(list('theta'=log(10),
                  'covariateProbs'=list(c(0.8,0.2),
                        c(0.2,0.8),
                        c(0.6,0.2,0.2),
                        c(0.25,0.5,0.25),
                        rep(0.25,4))),
            list('theta'=log(5),
                  'covariateProbs'=list(c(0.8,0.2),
                        c(0.9,0.1),
                        c(0.7,0.15,0.15),
                        c(0.6,0.2,0.2),
                        rep(0.25,4))),
            list('theta'=log(2.5),
                  'covariateProbs'=list(c(0.2,0.8),
                        c(0.9,0.1),
                        c(0.8,0.1,0.1),
                        c(0.8,0.1,0.1),
                        rep(0.25,4))),
            list('theta'=log(1),
                  'covariateProbs'=list(c(0.2,0.8),
                        c(0.3,0.7),
                        c(0.15,0.7,0.15),
                        c(0.8,0.1,0.1),
                        rep(0.25,4))),
            list('theta'=log(0.1),
                  'covariateProbs'=list(c(0.2,0.8),
                        c(0.5,0.5),
                        c(0.1,0.1,0.8),
                        c(0.8,0.1,0.1),
                        c(0.1,0.1,0.1,0.7)))))

clusSummaryPoissonNormal<-list(
   	'outcomeType'='Poisson',
   	'covariateType'='Normal',
   	'nCovariates'=2,
   	'nFixedEffects'=2,
		'fixedEffectsCoeffs'=c(-0.05,0.1),
   	'offsetLims'=c(0.9,1.1),
   	'missingDataProb'=0.001,
   	'nClusters'=3,
		'clusterSizes'=c(500,600,400),		
   	'clusterData'=list(list('theta'=log(10),
                  'covariateMeans'=c(0,2),
                  'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),                        
            list('theta'=log(3),
                  'covariateMeans'=c(3,2),
                  'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
            list('theta'=log(0.1),
                  'covariateMeans'=c(10,-5),
                  'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))            

clusSummaryBinomialNormal<-list(
   	'outcomeType'='Binomial',
   	'covariateType'='Normal',
   	'nCovariates'=2,
   	'nFixedEffects'=1,
		'fixedEffectsCoeffs'=c(0.1),
   	'nTrialsLims'=c(5,15),
   	'missingDataProb'=0.001,
   	'nClusters'=3,
		'clusterSizes'=c(200,700,100),		
   	'clusterData'=list(list('theta'=log(10),
                  'covariateMeans'=c(0,2),
                  'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),                        
            list('theta'=log(3),
                  'covariateMeans'=c(3,2),
                  'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
            list('theta'=log(0.1),
                  'covariateMeans'=c(10,-5),
                  'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))            

clusSummaryNormalNormal<-list(
   	'outcomeType'='Normal',
   	'covariateType'='Normal',
   	'nCovariates'=2,
   	'nFixedEffects'=1,
		'fixedEffectsCoeffs'=c(0.1),
   	'sigmaSqY'=1,
   	'missingDataProb'=0.001,
   	'nClusters'=3,
		'clusterSizes'=c(300,300,300),				
   	'clusterData'=list(list('theta'=-5,
                  'covariateMeans'=c(0,2),
                  'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),                        
            list('theta'=0,
                  'covariateMeans'=c(3,2),
                  'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
            list('theta'=5,
                  'covariateMeans'=c(10,-5),
                  'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))            


clusSummaryVarSelectBernoulliDiscrete<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=10,
      'nCategories'=rep(2,10),
      'nFixedEffects'=0,
      'missingDataProb'=0,
      'nClusters'=5,
		'clusterSizes'=c(200,200,200,200,200),				
      'clusterData'=list(list('theta'=log(1.0/9.0),
                  'covariateProbs'=list(c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9))),
            list('theta'=log(3.0/7.0),
                  'covariateProbs'=list(c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.1,0.9))),
            list('theta'=0,
                  'covariateProbs'=list(c(0.9,0.1),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9))),
            list('theta'=log(7.0/3.0),
                  'covariateProbs'=list(c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.1,0.9))),
            list('theta'=log(9.0),
                  'covariateProbs'=list(c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.1,0.9),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.9,0.1),
                        c(0.1,0.9)))))
