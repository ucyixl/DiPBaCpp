# Sample dataset
clusSummaryVS1<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=10,
      'isOrdinal'=rep(F,10),
      'nCategories'=rep(2,10),
      'nConfounders'=0,
      'missingDataProb'=0,
      'nClusters'=5,
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

clusSummaryVS2<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=7,
      'isOrdinal'=rep(F,7),
      'nCategories'=rep(3,7),
      'nConfounders'=0,
      'missingDataProb'=0,
      'nClusters'=3,
      'clusterData'=list(list('theta'=log(1.0/9.0),
                  'covariateProbs'=list(c(0.7,0.2,0.1),
                        c(0.7,0.2,0.1),
                        c(0.7,0.2,0.1),
                        c(0.7,0.2,0.1),
                        c(0.7,0.2,0.1),
                        c(0.7,0.2,0.1),
                        c(0.1,0.2,0.7))),
            list('theta'=0,
                  'covariateProbs'=list(c(0.1,0.2,0.7),
                        c(0.7,0.2,0.1),
                        c(0.1,0.2,0.7),
                        c(0.7,0.2,0.1),
                        c(0.1,0.2,0.7),
                        c(0.7,0.2,0.1),
                        c(0.1,0.2,0.7))),
            list('theta'=log(9.0),
                  'covariateProbs'=list(c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7),
                        c(0.1,0.2,0.7)))))


clusSummaryVS3<-list(
            'outcomeType'='Bernoulli',
            'covariateType'='Discrete',
            'nCovariates'=200,
            'isOrdinal'=rep(F,200),
            'nCategories'=rep(3,200),
            'nConfounders'=4,
            'confounderCoeff'=c(0.5,-0.5,0.3,1),
            'missingDataProb'=0,
            'nClusters'=5,
            'clusterData'=list(list('theta'=-2,
                        'covariateProbs'=vector(length=200,mode="list")),
                  list('theta'=-0.7,
                        'covariateProbs'=vector(length=200,mode="list")),
                  list('theta'=0,
                        'covariateProbs'=vector(length=200,mode="list")),
                  list('theta'=0.7,
                        'covariateProbs'=vector(length=200,mode="list")),
                  list('theta'=2,
                        'covariateProbs'=vector(length=200,mode="list"))))

p1<-c(0.9025,0.095,0.025)
p2<-c(0.3025,0.4950,0.2025)
p3<-c(0.5625,0.375,0.0625)
for(i in 1:200){
   if(i<=10){
      clusSummaryVS3$clusterData[[1]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[2]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[3]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[4]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[5]]$covariateProbs[[i]]<-p2
   }else if(i<=20){
      clusSummaryVS3$clusterData[[1]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[2]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[3]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[4]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[5]]$covariateProbs[[i]]<-p2      
   }else if(i<=30){
      clusSummaryVS3$clusterData[[1]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[2]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[3]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[4]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[5]]$covariateProbs[[i]]<-p2
   }else if(i<=40){
      clusSummaryVS3$clusterData[[1]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[2]]$covariateProbs[[i]]<-p2
      clusSummaryVS3$clusterData[[3]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[4]]$covariateProbs[[i]]<-p1
      clusSummaryVS3$clusterData[[5]]$covariateProbs[[i]]<-p2      
   }else{
      clusSummaryVS3$clusterData[[1]]$covariateProbs[[i]]<-p3
      clusSummaryVS3$clusterData[[2]]$covariateProbs[[i]]<-p3
      clusSummaryVS3$clusterData[[3]]$covariateProbs[[i]]<-p3
      clusSummaryVS3$clusterData[[4]]$covariateProbs[[i]]<-p3
      clusSummaryVS3$clusterData[[5]]$covariateProbs[[i]]<-p3      
   }
   
}      
      
      
clusSummary0<-list(
      'outcomeType'='Bernoulli',
      'covariateType'='Discrete',
      'nCovariates'=5,
      'isOrdinal'=rep(F,5),
      'nCategories'=c(3,3,3,3,3),
      'nConfounders'=0,
      'missingDataProb'=0.001,
      'nClusters'=5,
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



clusSummary1<-list(
      'outcomeType'='Poisson',
      'covariateType'='Discrete',
      'nCovariates'=5,
      'isOrdinal'=rep(F,5),
      'nCategories'=c(2,2,3,3,4),
      'nConfounders'=1,
      'confounderLinks'=c(1),
      'confounderStdDev'=c(0.25),
      'confounderCoeff'=c(0.01),
      'offsetLims'=c(0.9,1.1),
      'missingDataProb'=0.001,
      'nClusters'=5,
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

clusSummary2<-list(
   'outcomeType'='Poisson',
   'covariateType'='Normal',
   'nCovariates'=2,
   'nConfounders'=0,
   'offsetLims'=c(0.9,1.1),
   'missingDataProb'=0.001,
   'nClusters'=3,
   'clusterData'=list(list('theta'=log(10),
                              'covariateMeans'=c(0,2),
                              'covariateCovariance'=matrix(c(0.5,0,0,3),nrow=2)),                        
                        list('theta'=log(3),
                              'covariateMeans'=c(3,2),
                              'covariateCovariance'=matrix(c(1,0,0,1),nrow=2)),
                        list('theta'=log(0.1),
                              'covariateMeans'=c(10,-5),
                              'covariateCovariance'=matrix(c(2,0.7,0.7,1),nrow=2))))            
      
# Generate simulated data for testing C++ DiPBaC
generateXData<-function(nSubjects,subjectsPerCluster,clusterSummary,generateConfounders=F){

   nCovariates<-clusterSummary$nCovariates
   covariateType<-clusterSummary$covariateType
   missingDataProb<-clusterSummary$missingDataProb
   X<-matrix(-999,nSubjects,nCovariates)
   if(generateConfounders){
      nConfounders<-clusterSummary$nConfounders
      if(nConfounders>0){
         W<-matrix(0,nSubjects,nConfounders)
      }else{
         W<-NULL
      }
   }else{
      W<-NULL
   }
   
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
         isOrdinal<-clusterSummary$isOrdinal
         for(j in 1:nCovariates){
            if(i>1&&runif(1)<missingDataProb){
               X[i,j]<--999
            }else{
               if(!isOrdinal[j]){
                  u<-runif(1)
                  nCategories<-length(clusterData$covariateProbs[[j]])
                  for(kk in 1:nCategories){
                     if(u<cumsum(clusterData$covariateProbs[[j]])[kk]){
                        X[i,j]<-kk-1
                        break
                     }
                  }
               }else{
                  # This flips the ordinal model on its head to generate
                  # data, might also be a possibility for a revised ordinal
                  # model that would impose smoothing of probabilities
                  z<-rnorm(1,clusterData$covariateMeans[j],1)
                  varSet<-F
                  for(kk in 1:(nCategories-1)){
                     if(z<qnorm(kk/nCategories,0,1)){
                        X[i,j]<-kk-1
                        varSet<-T
                        break
                     }
                  }
                  if(!varSet){
                     # Must be in the last category
                     X[i,j]<-nCategories-1
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
      
      if(generateConfounders){
         # Loop over confounders
         if(nConfounders>0){
            for(j in 1:nConfounders){
               linkedCovariate<-clusterSummary$confounderLinks[j]
               stdDev<-clusterSummary$confounderStdDev[j]
               if(X[i,linkedCovariate]>=0){
                  W[i,j]<-rnorm(1,X[i,linkedCovariate],stdDev)
               }else{
                  # Just choose one of the past values and perturb it
                  # This can't happen for i=1 because missing values prevented
                  # for this individual when generating X
                  W[i,j]<-rnorm(1,W[as.integer(runif(1,1,i)),j],stdDev)
               }
            }
         }
      }
   }
   return(list('X'=X,'W'=W))
      
}

# Generate simulated data for testing C++ DiPBaC
generateYData<-function(nSubjects,subjectsPerCluster,W,clusterSummary){
   
   Y<-rep(0,nSubjects)
   outcomeType<-clusterSummary$outcomeType
   nConfounders<-clusterSummary$nConfounders
   if(nConfounders>0){
      beta<-clusterSummary$confounderCoeff   
   }
   if(outcomeType=='Poisson'){
      offset<-runif(nSubjects,clusterSummary$offsetLims[1],clusterSummary$offsetLims[2])
   }else{
      offset<-NULL
   }
   
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
      if(nConfounders>0){
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
      }
      
   }
   
   return(list("Y"=Y,"offset"=offset))
}

# Write the sample data to file
writeSampleDataToFile<-function(nSubjects,X,W,Y,offset,clusterSummary,fileName){

   write(nSubjects,fileName,append=F)
   nCovariates<-clusterSummary$nCovariates
   nConfounders<-clusterSummary$nConfounders
   write(nCovariates,fileName,append=T)
   covNames<-paste('Variable',seq(1,nCovariates,1),sep="")
   write(covNames,fileName,append=T,ncolumns=1)
   write(nConfounders,fileName,append=T)
   if(nConfounders>0){
      confNames<-paste('Confounder',seq(1,nConfounders,1),sep="")
      write(confNames,fileName,append=T,ncolumns=1)
   }
   if(clusterSummary$covariateType=='Discrete'){
      write(clusterSummary$nCategories,fileName,ncolumns=nCovariates,append=T)
      write(ifelse(clusterSummary$isOrdinal,1,0),fileName,ncolumns=nCovariates,append=T)
   }
   outData<-paste(Y,apply(X,1,paste,collapse=" "))
   if(nConfounders>0){
      outData<-paste(outData,apply(W,1,paste,collapse=" "))
   }
   if(clusterSummary$outcomeType=="Poisson"){
      outData<-paste(outData,offset)
   }
   write(outData,fileName,append=T,ncolumns=1)
}
