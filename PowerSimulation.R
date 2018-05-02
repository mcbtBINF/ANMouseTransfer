#Adapted from A. Fodor's code available on GitHub at : https://github.com/afodor/metagenomicsTools/blob/master/src/evanSimulations/simulatedPower.txt

# Power Simulation
rm(list=ls())

# this simulation assumes normality of data
numSimulations <- 100

# Number of taxa from phylum level to genus
numTaxa <- c(20, 36, 52, 129, 411)
numTaxa <- 411

numCases <- 40 #120 for mice
numControls <- 15 #40 in the case of paired T1vT2 in humans (120 for mice), 45 for HCvsAN 

# This value is calcuated from the number of significant hits surviving multiple hypothesis correction at the genus level for the mice at timepoint 4
fractionTruePositives <-c(0.05, 0.10, 0.20)
fractionTruePositives <- 0.05

# Standard lists of effect sizes
# In terms of Cohen's d: 0.10 = quite small, 0.2 = small, 0.50 = moderate, 0.80 = large
globalEffectSize <- c(0.10, 0.20, 0.50, 0.80)

globalEffectSize <- 0.24 # for T1 vs T2
globalEffectSize <- 0.43 # for AN vs HC

# The calculations connected to the background AN rate can be ignored, the "casecontrol" related calculation is the relevant one for this study
# https://www.mirror-mirror.org/eating-disorders-statistics.htm
# 3.7% with women and going to assume a similar, but slightly lower rate in men of 3.3%
fractionAN <- 0.037

casesAreAN <- rnorm(numCases) <= fractionAN 
controlsAreAN <- rnorm(numControls) <= fractionAN 
allAreAN <- c(casesAreAN, controlsAreAN)
allIsCase <- c( rep(TRUE, numCases), rep(FALSE, numControls) )

powerCaseControl <- vector()
powerAN <- vector()
powerInteraction <- vector()

for( i in 1:numSimulations)
{
  pValuesCaseControl <- vector()
  pValuesAN <- vector()
  pValuesInteraction <- vector()
  taxaisCaseControlPositive <- rnorm(numTaxa) <= fractionTruePositives
  taxaisANPositive <- rnorm(numTaxa) <= fractionTruePositives
  
  for( j in 1:numTaxa)
  {
    controlVals <- vector()
    
    for ( k in 1:numControls ) 
    {
      thisEffectSize = 0;
      if( controlsAreAN[k] &  taxaisANPositive[j] ) 
      {
        thisEffectSize = globalEffectSize;
      }
      
      controlVals[k] = rnorm(1, mean= thisEffectSize )
    }
    
    caseVals <- vector()
    
    for ( k in 1:numCases ) 
    {
      thisEffectSize = 0;
      
      if( taxaisCaseControlPositive[j] ) 
      {
        thisEffectSize = globalEffectSize;
      }
      
      if( casesAreAN[k] &  taxaisANPositive[j] ) 
      {
        thisEffectSize = thisEffectSize +  globalEffectSize;
      }
      
      if( casesAreAN[k] &  taxaisANPositive[j] &  taxaisCaseControlPositive[j] ) 
      {
        # if the taxa is responsive to both case and control and AN
        # there is an interaction doubling the impact of either one alone
        thisEffectSize = thisEffectSize + 2 * globalEffectSize
      }
      
      caseVals[k] = rnorm(1, mean= thisEffectSize )
    }
    
    # horribly inefficient here, but seems to still be fast enough
    data <- c( caseVals, controlVals ) 
    
    # build model with the interaction term
    myLm <- lm( data ~ allIsCase * allAreAN ) 
    
    #record p-values with anova
    myAnova <- anova(myLm)
    
    # use ANOVA to capture the p-values
    pValuesCaseControl[j] <- anova(myLm)$"Pr(>F)"[1]
    pValuesAN[j] <- anova(myLm)$"Pr(>F)"[2]
    pValuesInteraction[j] <- anova(myLm)$"Pr(>F)"[3]
  }
  
  # adjust the p-values
  caseControlAdjust <- p.adjust( pValuesCaseControl, method="BH")
  
  numTrueFound <-0 
  
  for ( j in 1:numTaxa ) 
  {
    if( caseControlAdjust[j] < 0.10 & taxaisCaseControlPositive[j] )
      numTrueFound = numTrueFound + 1
  }
  
  powerCaseControl[i] =  numTrueFound  / sum(taxaisCaseControlPositive)
  
  ANAdjust <- p.adjust( pValuesAN, method="BH")
  
  numTrueFound <-0 
  
  for ( j in 1:numTaxa ) 
  {
    if( ANAdjust[j] < 0.10 & taxaisANPositive[j] )
      numTrueFound = numTrueFound + 1
  }
  
  powerAN[i] =  numTrueFound  / sum(taxaisANPositive)
  
  interactionAdjust <- p.adjust( pValuesInteraction, method="BH")
  
  numTrueFound <-0 
  
  for ( j in 1:numTaxa ) 
  {
    if( interactionAdjust[j] < 0.10 & taxaisANPositive[j] &  taxaisCaseControlPositive[j] )
      numTrueFound = numTrueFound + 1
  }
  
  powerInteraction[i] =  numTrueFound  / sum(taxaisANPositive &taxaisCaseControlPositive )
  
  
}

# print out mean power across all simulations
print(c(numCases, numControls, globalEffectSize, fractionTruePositives))
print(c(mean(powerCaseControl),mean(powerAN),mean(powerInteraction)))
