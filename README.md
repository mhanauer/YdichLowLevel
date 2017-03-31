# YdichLowLevel
# Jags-Ydich-XmetSsubj-Dicot-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
# This has the most updated version but won't work
source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , x2Name = "x2",x3Name = "x3", yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  x2 = data[,x2Name]
  x3 = data[,x3Name]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    x2 = x2,
    x3 = x3,
    y = y ,
    s = s ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
  Ntotal <- length(y)
  xm <- mean(x)
  x2m <- mean(x2)
  x3m <- mean(x3)
  xsd <- sd(x)
  x2sd <- sd(x2)
  x3sd <- sd(x3)
  
  for ( i in 1:length(y) ) {
  zx[i] <- ( x[i] - xm ) / xsd
  zx2[i] <- ( x2[i] - x2m ) / x2sd
  zx3[i] <- ( x3[i] - x3m ) / x3sd
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dbern(mu[i]) 
  mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx2[i] + zbeta3[s[i]] * zx3[i]))
  }
  guess ~ dbeta(1,9)
  for ( j in 1:Nsubj ) {
  zbeta0[j] ~ dnorm( 0 , 1/(10)^2  )  
  zbeta1[j] ~ dnorm( 0 , 1/(10)^2  )
  zbeta2[j] ~ dnorm( 0 , 1/(10)^2  )
  zbeta3[j] ~ dnorm( 0 , 1/(10)^2  )


  }
  # Priors vague on standardized scale:
  zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta2mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta3mu ~ dnorm( 0 , 1/(10)^2 )
  
  # Transform to original scale:
  for ( j in 1:Nsubj ) {
  beta1[j] <- zbeta1[j]  / xsd 
  beta2[j] <- zbeta2[j]  / x2sd
  beta3[j] <- zbeta3[j]  / x3sd
  beta0[j] <- zbeta0[j]  - zbeta1[j] * xm / xsd + zbeta2[j] * x2m / x2sd + zbeta3[j] * x3m / x3sd
  }
  beta1mu <- zbeta1mu  / xsd
  beta2mu <- zbeta2mu  / x2sd
  beta3mu <- zbeta3mu  / x3sd
  beta0mu <- zbeta0mu - zbeta1mu * xm / xsd + zbeta2mu * x2m / x2sd + zbeta3mu * x3m  / x3sd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1","beta2" , "beta3", "beta0mu" , "beta1mu" , "beta2mu", "beta3mu",
                  "zbeta0" , "zbeta1" , "zbeta2", "zbeta3",  "zbeta0mu" , "zbeta1mu" ,"zbeta2mu", "zbeta3mu")
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}
