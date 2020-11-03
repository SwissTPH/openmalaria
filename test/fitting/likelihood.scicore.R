# Revision of OpenMalaria_calculateLF_2016.R
#
# Takes output from 61 scenarios (used in initial fitting to 11 objectives for OpenMalaria) and
# 1) Calculates loss functions for each fitting objective (RSS, loglikelihoods)
# 2) Plots predictions against observed data
# 
# Authors
# Melissa Penny (May 2020): Fix errors introduced by expected severe and deaths and MOI 
# Ewan Cameron (2016): minor i/o mods & pretty plots
# Melissa Penny (2012-2016): overall structure, data, obtaining predictions, output, naming conventions, plots, fix errors from phase A, calculate and use log-likelihoods (and not RSS from 2006 and Plos Med paper) 
# Thomas Smith (2016): updated for v35 expected severe and deaths
# Initial plotting from PhaseA(2006) by AR: 6.7.07
# New file (with errors fixed) MP: 14.05.2012 

library(gridExtra)

# =============================================================================
# CODES/IDS FIELD DATA AND MODEL RESULTS codes/ids for for field data and model predictions
# =============================================================================
# Function fieldData_nameCodes -  In the field data, the codes used for the different measures are as follows:
fieldData_nameCodes <- function(){
	fieldDataCode <- list()
    fieldDataCode$nHost <- 0;  #number of hosts
    fieldDataCode$nInfect <- 1;  #number of infected hosts
    fieldDataCode$nExpectd <- 2;  #expected number of infected hosts
    fieldDataCode$npatent <- 3;  #number of patent hosts
    fieldDataCode$sumlogdens <- 5;  #sum of the logarithm of the density
    fieldDataCode$pyr <- 10;  #person-years
    fieldDataCode$totalpatentInf <- 13;  #total patent infections
   	fieldDataCode$malMort <- 16;  #Malaria mortality rate
    fieldDataCode$recInc <- 18;  #recorded incidence disease
    fieldDataCode$pyrogenThres <- 19;  #pyrogenic threshold as para/leuc ratio
    fieldDataCode$rate <- 20;  #relative risk of severe compared to agegp 1
    fieldDataCode$prop <- 21;  #proportion (unspecified)
    fieldDataCode$allCauseIMR <- 22;  #all-cause infant mortality rate
    return(fieldDataCode)
}

# Function modelResults_nameCodes -  for model results from openMalaria, the codes used for the different measures are as follows:
modelResults_nameCodes <- function(){
	modelResultsCode <- list()		
	modelResultsCode$nHost <- 0;   #Total number of humans
	modelResultsCode$nInfect <- 1;   #number of infected hosts
	modelResultsCode$nExpectd <- 2;   #expected number of infected hosts
	modelResultsCode$nPatent <- 3;   #number of patent hosts
	modelResultsCode$sumLogPyrogenThres <- 4;   #Sum of the log of the pyrogen threshold
	modelResultsCode$sumlogDens <- 5;   #Sum of the logarithm of the parasite density
	modelResultsCode$totalInfs <- 6;   #Total infections
	modelResultsCode$nTransmit <- 7;   # Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes. Single value, not per age-group.
	modelResultsCode$totalPatentInf <- 8;   #Total patent infections
	modelResultsCode$sumPyrogenThresh <- 10;   #Sum of the pyrogenic threshold
	modelResultsCode$nTreatments1 <- 11;   #number of treatments (1st line) (added to 1-day model in 24.1)
	modelResultsCode$nTreatments2 <- 12;   #number of treatments (2nd line) (added to 1-day model in 24.1)
	modelResultsCode$nTreatments3 <- 13;   #number of treatments (inpatient) (added to 1-day model in 24.1)
	modelResultsCode$nUncomp <- 14;   #number of episodes (uncomplicated)
	modelResultsCode$nSevere <- 15;   #number of episodes (severe)
	modelResultsCode$nSeq <- 16;   #recovered cases with sequelae
	modelResultsCode$nHospitalDeaths <- 17;   #deaths in hospital
	modelResultsCode$nIndDeaths <- 18;   #number of deaths (indirect)
	modelResultsCode$nDirDeaths <- 19;   #number of deaths (direct)
	modelResultsCode$nEPIVaccinations <- 20;   #number of EPI vaccine doses given
	modelResultsCode$allCauseIMR <- 21;   #all cause infant mortality rate (returned as a single number over whole intervention period, instead of from a survey interval)
	modelResultsCode$nMassVaccinations <- 22;   #number of Mass / Campaign vaccine doses given
	modelResultsCode$nHospitalRecovs <- 23;   #recoveries in hospital without sequelae
	modelResultsCode$nHospitalSeqs <- 24;   #recoveries in hospital with sequelae
	modelResultsCode$nIPTDoses <- 25;   #number of IPT Doses
	modelResultsCode$annAvgK <- 26;   #Annual Average Kappa. Calculated once a year as sum of human infectiousness weighted by initial EIR for that time of year.
	modelResultsCode$nNMFever <- 27;   #Number of episodes of non-malaria fever
	modelResultsCode$expectedDirectDeaths <- 74;   #Expected number of deaths (direct)
	modelResultsCode$expectedHospitalDeaths <- 75;   #Expected number of deaths in hospital
	modelResultsCode$expectedIndirectDeaths <- 76;   #Expected number of deaths (indirect)
	modelResultsCode$expectedSevere <- 78;   #Expected number of episodes (severe)
 	return(modelResultsCode)
}

# =============================================================================
# FIELD DATA
# =============================================================================
# Incidence Field Data (scenarios 30 Matsari after intervention)
incidenceAfterInterventionFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Matsari")
	field$site <- site
	field$scenarioNum <- c(30)
	age<-c(0.25,0.75,1.5,3,5,8,12.5,17.5,25,35,50,80) # age
	for (i in 1:1) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	# get data and calculate prevalence and density
	for (i in 1:1) {
		scNum <- field$scenarioNum[i]
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){	
			# data over all surveys
			field$data[[site[i]]]$nPatent[ageGroup] <- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$npatent & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$nHosts[ageGroup] <- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$prevalence[ageGroup] <- field$data[[site[i]]]$nPatent[ageGroup]/field$data[[site[i]]]$nHosts[ageGroup]
			# data for individual surveys
			field$data[[site[i]]]$indivSurveys_nPatent[[ageGroup]]<- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$npatent & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]] <- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_prevalence[[ageGroup]] <- field$data[[site[i]]]$indivSurveys_nPatent[[ageGroup]]/field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]]
			
		}
	}
					
	return(field)
}
	
# Asexual Field Data (scenarios 24,28,29,35,34 and 31)
asexualFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Sugungum","Rafin Marke","Matsari","Namawala","Navrongo","Idete")
	field$site <- site
	field$scenarioNum <- c(24,28,29,35,34,31)
	age<-c(0.25,0.75,1.5,3,5,8,12.5,17.5,25,35,50,80) # age
	for (i in 1:5) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	field$data[[site[6]]]$age <- c(1/12, 0.25, 0.416666667, 0.583333333, 0.75, 0.916666667, 1.5, 2.5, 3.5, 4.5, 5.5)
	field$data[[site[6]]]$name <- site[6]
	
	# get data and calculate prevalence and density
	for (i in 1:6) {
		scNum <- field$scenarioNum[i]
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){	
			# data over all surveys
			field$data[[site[i]]]$nHosts[ageGroup] <- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$nPatent[ageGroup] <- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$npatent & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$sumlogDens[ageGroup]<- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$sumlogdens & fieldDataAll$V3== ageGroup])
			
			field$data[[site[i]]]$prevalence[ageGroup] <- field$data[[site[i]]]$nPatent[ageGroup]/field$data[[site[i]]]$nHosts[ageGroup]
			field$data[[site[i]]]$density[ageGroup] <- exp(field$data[[site[i]]]$sumlogDens[ageGroup]/field$data[[site[i]]]$nPatent[ageGroup])
			
			# data for individual surveys
			field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]] <- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_nPatent[[ageGroup]]<- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$npatent & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_sumlogDens[[ageGroup]]<- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$sumlogdens & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_prevalence[[ageGroup]] <- field$data[[site[i]]]$indivSurveys_nPatent[[ageGroup]]/field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]]
			field$data[[site[i]]]$indivSurveys_density[[ageGroup]] <- exp(field$data[[site[i]]]$indivSurveys_sumlogDens[[ageGroup]]/field$data[[site[i]]]$indivSurveys_nPatent[[ageGroup]])
		}
	}
	return(field)
}	

# moi Field Data (scenario 34)
moiFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Navrongo")
	field$site <- site
	field$scenarioNum <- c(34)
	field$data[[site[1]]]$age <- c(0.25,0.75,1.5,3,5,8,12.5,17.5,25,35,50,80) # age

	for (i in 1:1) {
		scNum <- field$scenarioNum[i]
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){
		# data over all surveys
			field$data[[site[i]]]$nHosts[ageGroup]<- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$numPatInf[ageGroup]<- sum(fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$totalpatentInf & fieldDataAll$V3== ageGroup])
			field$data[[site[i]]]$moi[ageGroup]<- field$data[[site[i]]]$numPatInf[ageGroup]/field$data[[site[i]]]$nHosts[ageGroup]
		# data for individual surveys
			field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]] <- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$nHost & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_numPatInf[[ageGroup]]<- fieldDataAll$V5[fieldDataAll$V1==scNum &fieldDataAll$V4==fieldDataID$totalpatentInf & fieldDataAll$V3== ageGroup]
			field$data[[site[i]]]$indivSurveys_moi[[ageGroup]] <- field$data[[site[i]]]$indivSurveys_numPatInf[[ageGroup]]/field$data[[site[i]]]$indivSurveys_nHosts[[ageGroup]]
		}		
	}	
	return(field)
}	
		
# acute Field Data (scenarios 232 and 233)
acuteSenegalFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Ndiop","Dielmo")
	field$site <- site
	field$scenarioNum <- c(232, 233)
	age<-c(seq(1:15)+0.5,17.5,22.5,27.5,35,45,55,80)
	for (i in 1:2) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	# only use survey 2
	surNum <- 2# outcome 18 is the recorded incidence 
	for (i in 1:2) {
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){	
			scNum <- field$scenarioNum[i]
			field$data[[site[i]]]$incid[ageGroup]<- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V2==surNum & fieldDataAll$V4==fieldDataID$recInc & fieldDataAll$V3== ageGroup]
		}
	}	
	return(field)
}

# acute Idete Field Data (scenario 49)
acuteIdeteFieldData <- function(fieldDataAll, fieldDataID){	
	field <- list()
	site<-c("Idete")
	field$site <- site
	field$scenarioNum <- c(49)
	age <- c(0.125,0.375,0.625,0.875)
	for (i in 1:1) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	# outcome 18 is the recorded incidience 
	for (i in 1:1) {
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){	
			scNum <- field$scenarioNum[i]
			field$data[[site[i]]]$incid[ageGroup]<- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V4==fieldDataID$recInc & fieldDataAll$V3== ageGroup]
		}
	}	
	return(field)
}

# Acute Idete Field Data (scenario 234)
acutePyGenThresFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("NdiopDielmo")
	field$site <- site
	field$scenarioNum <- c(234)
	age <-c(seq(1:15)+0.5,17.5,22.5,27.5,35,45,55,80)
	for (i in 1:1) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	# outcome 19 is the recorded incidience 
	for (i in 1:1) {
		#field$data[[site[i]]]$Ystar[1]<-NA # no data for age group 1
		for (ageGroup in 1:length(field$data[[site[i]]]$age)){	
			scNum <- field$scenarioNum[i]
			field$data[[site[i]]]$Ystar[ageGroup]<- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V4==fieldDataID$pyrogenThres & fieldDataAll$V3== ageGroup]
		}
		# no data for age group set Ystart to NA
		field$data[[site[i]]]$Ystar[field$data[[site[i]]]$Ystar<0] <- NA   
	}	
	return(field)
}

# Severe Field Data (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)
severeMarshSnowFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	field$scenarioNum <- c(501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527) # scenario numbers
	# Marsh-Snow data - the field data is the same for all the scenarios, so just take scenario 501
	scNum <- field$scenarioNum[1]
	# access to treatment fitted from phase A
	accessSevereTreatment <- 0.48
	field$MarshSnow_episodes <- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V4==fieldDataID$rate ]/accessSevereTreatment
	field$prevalence <- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V4==fieldDataID$prop ]*100 # for plotting
	
	return(field)
}

# Age spcific prevalance - severe data (scenario 158,167,173,176)
severeSnow1997FieldData <- function(fieldDataAll, fieldDataID){	
	field <- list()
	site<-c("Sukuta","Kilifi-north","Kilifi-south","Siaya")
	field$site <- site
	field$scenarioNum <- c(158,167,173,176) # scenario numbers
	# data from Snow et al 1997, table 2 for 4 of 5 sites
	age <- c(0.5,3,7) # age
	for (i in 1:length(field$scenarioNum)) {
		field$data[[site[i]]]$age <- age
		field$data[[site[i]]]$name <- site[i]
	}
	# access to treatment fitted from phase A
	accessSevereTreatment <- 0.48
	# get data 
	for (i in 1:length(field$scenarioNum)) {
		for (ageGroup in 2:length(field$data[[site[i]]]$age)){	
			scNum <- field$scenarioNum[i]
			field$data[[site[i]]]$relativeRiskComparedAgeGroup1[ageGroup]<- sum(fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V3== ageGroup &fieldDataAll$V4==fieldDataID$rate])
			
		}
	}
	# data from Snow et al 1997, table 2 for 4 sites
	# original data for plotting, this is not in the fielddata.txt file
	# values are (1-11 months, 1-4 years, 5-9 years)
	field$data[[site[1]]]$rate_epPer1000perYr <-c(23.3,35.3,16.3)/accessSevereTreatment
	field$data[[site[2]]]$rate_epPer1000perYr <-c(59.5,41.7,5.3)/accessSevereTreatment
	field$data[[site[3]]]$rate_epPer1000perYr <-c(79.0,17.4,1.7)/accessSevereTreatment
	field$data[[site[4]]]$rate_epPer1000perYr <-c(84.6,18.8,1.7)/accessSevereTreatment
	return(field)
}

# Direct mortality (scenarios 301,302,303,312,316,317,318,326,327) 
directMortalityFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Bo, Sierra Leone","Niakhar, Senegal","Farafenni South Bank, The Gambia","Kilifi north, Kenya","Navrongo, Ghana","Asembo Bay, Kenya","Bagamoyo/Yombo Tanzania","Bandafassi, Senegal","Kongondjan area , Burkina Faso17")
	field$site <- site
	field$scenarioNum <- c(301,302,303,312,316,317,318,326,327) # scenario numbers
	# data from korenromp 2003 table 1
	for (i in 1:length(field$scenarioNum)) {
		field$data[[site[i]]]$name <- site[i]
	}
	# get data 
	# data from korenromp 2003 table 1 - from table directly, not in text file used for fitting
	EIR <- c(34.606,11.675,8.9498,10.518,404.5,238.7,219.82,365.11,133.55)	
	#fieldEst <- c(12.8,10.9,9.4,9.2,9.3,20.8,22.1,4.7,2.2) # field estimates of under 5y mortality (these are values in text file multiplied by 1000 to get to rate per 1000 per year)
	fieldEstUpperCI <-c(21.2,15.6,14.6,12.9,17,29.8,33.1,7.5,5.2) # upper CI of field estimate	
	fieldEstLowerCI <-c(4.8,6.9,4.8,5.9,1.9,12.8,14,2.1,0) # lower CI	of field estimate

	for (i in 1:length(field$scenarioNum)) {
		ageGroup <-1
		scNum <- field$scenarioNum[i]
		field$data[[site[i]]]$dirDeathRateForAgeGroup1[ageGroup]<- sum(fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V3==ageGroup &fieldDataAll$V4==fieldDataID$malMort])
		field$data[[site[i]]]$EIR[ageGroup]<- EIR[i]
		field$data[[site[i]]]$dirDeathRate_epPer1000perYr[ageGroup]<- field$data[[site[i]]]$dirDeathRateForAgeGroup1[ageGroup]*1000
		field$data[[site[i]]]$dirDeathRate_epPer1000perYr_UpperCI[ageGroup]<- fieldEstUpperCI[i]
		field$data[[site[i]]]$dirDeathRate_epPer1000perYr_LowerCI[ageGroup]<- fieldEstLowerCI[i]
	}
	return(field)
}	

# Indirect mortality (scenarios 401,402,408,411,414,415,416,417,418,422,426) 
indirectMortalityFieldData <- function(fieldDataAll, fieldDataID){
	field <- list()
	site<-c("Bo, Sierra Leone", "Niakhar, Senegal", "Upper River Divison, The Gambia", "Karangasso, Burkina Faso", "Manhica, Mozambique", "Namawala, Tanzania", "Navrongo, Ghana", "Saradidi, Kenya", "Yombo, Tanzania", "Mlomp, Senegal", "Bandafassi, Senegal")
	
	# this data from table in paper isn't used anymore? Area I-V The Gambia (1992)
	field$site <- site
	field$scenarioNum <- c(401,402,408,411,414,415,416,417,418,422,426) # scenario numbers
	for (i in 1:length(field$scenarioNum)) {
		field$data[[site[i]]]$name <- site[i]
	}
	# get data 
	# data from many sources, summarised in Ross 2006
	EIR <- c(34.606,11.675,9.5048,244.15,38.402,329.08,404.5,238.7,219.82,30.176,365.11)
	for (i in 1:length(field$scenarioNum)) {
		ageGroup <-1
		scNum <- field$scenarioNum[i]
		field$data[[site[i]]]$indirDeath[ageGroup]<- fieldDataAll$V5[fieldDataAll$V1==scNum & fieldDataAll$V3==ageGroup &fieldDataAll$V4==fieldDataID$allCauseIMR]
		field$data[[site[i]]]$EIR[ageGroup]<- EIR[i]	
	}
	return(field)
}	

# =============================================================================
# EXTRACT MODEL PREDICTIONS - extract data from model predictions
# =============================================================================

# Function to allow choice of predictions between expected numbers of events or counts of discrete events
preferredOutput <- function(output, surNum, ageGroup, preferredMeasure, reserveMeasure){
  value <- NA
  reserve <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== reserveMeasure]
  preferred <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== preferredMeasure]
  if (length(reserve) > 0) value <- reserve 
  if (length(preferred) > 0) value <- preferred
  return(value)
}

# Function extractIncid - extract from output files prevalence and densites and sum over all surveys for each age group 
extractIncid<-function(scenNum, runDir, modelResultID) {
	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))
	
	# find output number of host, number of patent hosts and  Sum of the logarithm of the parasite density
	pred <- list()
	numAgeGroups <- max(output$V2)
	for (ageGroup in 1:numAgeGroups) {
		# predictions over all surveys
		pred$nHosts[ageGroup] <- sum(output$V4[output$V2==ageGroup & output$V3==modelResultID$nHost])
		pred$nPatent[ageGroup] <- sum(output$V4[output$V2==ageGroup  & output$V3==modelResultID$nPatent])
		pred$prevalence[ageGroup] <- pred$nPatent[ageGroup]/pred$nHosts[ageGroup]
		# predictions for individual surveys
		pred$indivSurveys_nHosts[[ageGroup]] <- output$V4[output$V2==ageGroup & output$V3==modelResultID$nHost]
		pred$indivSurveys_nPatent[[ageGroup]] <- output$V4[output$V2==ageGroup  & output$V3==modelResultID$nPatent]
		pred$indivSurveys_prevalence[[ageGroup]] <- pred$indivSurveys_nPatent[[ageGroup]]/pred$indivSurveys_nHosts[[ageGroup]]
	}	
	return(pred)
}

# Function extractAS - extract from output files prevalence and densites and sum over all surveys for each age group 
extractAS<-function(scenNum,runDir, densBias, modelResultID) {
	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))
	
	# find output number of host, number of patent hosts and  Sum of the logarithm of the parasite density
	pred <- list()
	numAgeGroups <- max(output$V2)
	for (ageGroup in 1:numAgeGroups) {
		# predictions over all surveys
		pred$nHosts[ageGroup] <- sum(output$V4[output$V2==ageGroup & output$V3==modelResultID$nHost])
		pred$nPatent[ageGroup] <- sum(output$V4[output$V2==ageGroup  & output$V3==modelResultID$nPatent])
		pred$sumlogDens[ageGroup] <- sum(output$V4[output$V2==ageGroup  & output$V3==modelResultID$sumlogDens])
		pred$prevalence[ageGroup] <- pred$nPatent[ageGroup]/pred$nHosts[ageGroup]
		pred$density[ageGroup] <- exp(pred$sumlogDens[ageGroup]/pred$nPatent[ageGroup])/densBias
		
		# predictions for individual surveys
		pred$indivSurveys_nHosts[[ageGroup]] <- output$V4[output$V2==ageGroup & output$V3==modelResultID$nHost]
		pred$indivSurveys_nPatent[[ageGroup]] <- output$V4[output$V2==ageGroup  & output$V3==modelResultID$nPatent]
		
		pred$indivSurveys_sumlogDens[[ageGroup]] <- output$V4[output$V2==ageGroup  & output$V3==modelResultID$sumlogDens]
		pred$indivSurveys_prevalence[[ageGroup]] <- pred$indivSurveys_nPatent[[ageGroup]]/pred$indivSurveys_nHosts[[ageGroup]]
		pred$indivSurveys_density[[ageGroup]] <- exp(pred$indivSurveys_sumlogDens[[ageGroup]]/pred$indivSurveys_nPatent[[ageGroup]])/densBias
		}	
	return(pred)
}

# Function extractMoi-extract from output files moi and sum over all surveys for each age group 
extractMoi <-function(scenNum,runDir, modelResultID){
	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))
	
	# find output required # 8 is totalPatentInf, 0 = nHosts
	pred <- list()
	numAgeGroups <- max(output$V2)
	for (ageGroup in 1:numAgeGroups) {
		# predictions over all surveys
		pred$nHosts[ageGroup] <- sum(output$V4[output$V2==ageGroup & output$V3== modelResultID$nHost])
		pred$numPatInf[ageGroup] <- sum(output$V4[output$V2==ageGroup  & output$V3== modelResultID$totalPatentInf])
		pred$moi[ageGroup] <- pred$numPatInf[ageGroup]/pred$nHosts[ageGroup]
		
		# predictions for individual surveys
		pred$indivSurveys_nHosts[[ageGroup]] <- output$V4[output$V2==ageGroup & output$V3==modelResultID$nHost]
		pred$indivSurveys_numPatInf[[ageGroup]] <- output$V4[output$V2==ageGroup  & output$V3==modelResultID$totalPatentInf]
		pred$indivSurveys_moi[[ageGroup]] <- pred$indivSurveys_numPatInf[[ageGroup]]/pred$indivSurveys_nHosts[[ageGroup]]
	}	
	return(pred)
}

# Function extractAEpSenegal - extract acute Episodes estimates 
extractAEpSenegal<-function(scenNum,runDir, bias, modelResultID) {	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))	
	# find output required #14 = uncomplicated, 15 = severe, 0 =nHosts
	pred <- list()
	numAgeGroups <- max(output$V2)
	# use only survey 2 (post intervention)
	surNum <- 2
	for (ageGroup in 1:numAgeGroups) {
		pred$nHosts[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nHost]
		pred$uncomp[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nUncomp]
		# if predictions for expected numbers of severe events are available use these.  Otherwise use the number of events
		pred$severe[ageGroup] <- preferredOutput(output, surNum, ageGroup, modelResultID$expectedSevere, modelResultID$nSevere)
    #pred$severe[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nSevere]
    pred$incid[ageGroup] <- ((pred$uncomp[ageGroup]+pred$severe[ageGroup])/pred$nHosts[ageGroup])/bias
	}	
	return(pred)
}
# this is for testing only TODO
extractAEpSenegal_survey1<-function(scenNum,runDir, bias, modelResultID) {	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))	
	# find output required #14 = uncomplicated, 15 = severe, 0 =nHosts
	pred <- list()
	numAgeGroups <- max(output$V2)
	# use only survey 2 (normally)
	# TODO check only survey 1
	surNum <- 1
	for (ageGroup in 1:numAgeGroups) {
		pred$nHosts[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nHost]
		pred$uncomp[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nUncomp]
		# if predictions for expected numbers of severe events are available use these.  Otherwise use the number of events
		pred$severe[ageGroup] <- preferredOutput(output, surNum, ageGroup, modelResultID$expectedSevere, modelResultID$nSevere)
		#pred$severe[ageGroup] <- output$V4[output$V1==surNum & output$V2==ageGroup & output$V3== modelResultID$nSevere]
		pred$incid[ageGroup] <- ((pred$uncomp[ageGroup]+pred$severe[ageGroup])/pred$nHosts[ageGroup])/bias
	}	
	return(pred)
}

# Function extractAEpIdete - extract acute Episodes estimates 
# extract incidence for Idete
extractAEpIdete <-function(scenNum,runDir, bias, modelResultID) {	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))	
	# find output required #14 = uncomplicted (not severe), 0 =nHosts
	pred <- list()
	numAgeGroups <- max(output$V2)
	for (ageGroup in 1:numAgeGroups) {
		pred$nHosts[ageGroup] <- output$V4[output$V2==ageGroup & output$V3== modelResultID$nHost]
		pred$uncomp[ageGroup] <- output$V4[output$V2==ageGroup & output$V3== modelResultID$nUncomp]
		pred$incid[ageGroup] <- ((pred$uncomp[ageGroup])/pred$nHosts[ageGroup])/bias

	}	
	return(pred)
}

# Function extractPyGen - extract pyrogenThres estimates 
extractPyGen <-function(scenNum,runDir, bias, modelResultID) {	
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_",scenNum,".txt", sep=""))	
	# find output required #4 = sumLogPyrogenThres, 0 =nHosts
	# we only use survey 1 as we want predictions from the same season as the field cases- but we have a high tolerance of what counts as the same season, i.e. the same quarter.
	surv <- 1
	pred <- list()
	numAgeGroups <- max(output$V2)
	for (ageGroup in 1:numAgeGroups) {
		pred$nhost[ageGroup]<-sum(output$V4[output$V1==surv & output$V2== ageGroup & output$V3== modelResultID$nHost])	
		pred$logpyrogt[ageGroup]<-sum(output$V4[output$V1==surv & output$V2== ageGroup & output$V3== modelResultID$sumLogPyrogenThres])
		pred$pyrogt[ageGroup]<-(exp(pred$logpyrogt[ageGroup]/pred$nhost[ageGroup]))*(bias)
		# which is exp(predlogP+log(bias))
	}
			
	return(pred)
}

# Function extractSM - extract SEVERE MALARIA estimates (assumes survey times to 24 months) 
extractSM<-function(scenNum,runDir,modelResultID) { 
	# read in predictions for particular scenario number
	outputTemp<-read.table(paste(runDir,"wu_5000_", scenNum,".txt", sep="")) 
	output<-outputTemp[outputTemp$V2==1,] # first age-group only (0-9yrs) 
	pred <- list()
	pred$nhost<-mean(output$V4[output$V3== modelResultID$nHost]) # mean num hosts
	pred$npatent<-mean(output$V4[output$V3== modelResultID$nPatent]) # mean num patent
	pred$prevalence<-(pred$npatent/pred$nhost)*100 
	# multiply by 100 to get to percentage
	# if predictions for expected numbers of severe events are available use these.  Otherwise use the number of events
	pred$severe<- 0
  	for (surNum in 1:24) {
    	pred$severe<- pred$severe + preferredOutput(output, surNum, 1, modelResultID$expectedSevere, modelResultID$nSevere)
  	}
	#pred$severe<-sum(output$V4[output$V3== modelResultID$nSevere]) 
	#require rate episodes per 1000 person years (so divide pred$severe/2)
	pred$sevEpPer1000personYr<-(pred$severe/(2*pred$nhost))*1000 
	return(pred) 
}
 
# Function extractSMRR -extract SEVERE MALARIA Relative Risk predictions (assuming single survey after 2 years)
extractSMRR<-function(scenNum,runDir, modelResultID) { 
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_", scenNum,".txt", sep="")) 
	# output is for 4 age groups (0-1, 1-5, 5-10, 10-99) at survey at 2 years
	pred <- list()
	pred$nhost<-output$V4[output$V3==modelResultID$nHost]
	# pred$severe <-output$V4[output$V3==modelResultID$nSevere]
	
	# if predictions for expected numbers of severe events are available use these.  Otherwise use the number of events
	pred$severe<- array(0, dim =c(24*4))
  	for (surNum in 1:24) {
  		for (AgeNum in 1:4) {
    		pred$severe[surNum*4-(4-AgeNum)]<- preferredOutput(output, surNum, AgeNum, modelResultID$expectedSevere, modelResultID$nSevere)
  		}
  	}
  	
	# require severe episodes per 1000 person years, so multiply by 1000 and divide by 2
	pred$sevEpPer1000personYr<-(pred$severe/(2*pred$nhost))*1000 
	return(pred)
}

extractMort <-function(scenNum,runDir, modelResultID) { 
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_", scenNum,".txt", sep="")) 
	# output is for 4 age groups (0-5, 5-10, 10-20, 20-99) at survey at 2 years, only require age group 1 (0-5 years)
	pred <- list()
	ageGroup <- 1
	pred$nhost<-output$V4[output$V2== ageGroup & output$V3==modelResultID$nHost]
	# if predictions for expected numbers of deaths are available use these.  Otherwise use the number of events
	pred$dirDeaths<- preferredOutput(output, 1, ageGroup, modelResultID$expectedDirectDeaths, modelResultID$nDirDeaths)
  #	pred$dirDeaths <-output$V4[output$V2== ageGroup &output$V3==modelResultID$nDirDeaths]
	# require direct deaths episodes per 1000 person years, so multiply by 1000 and divide by 2
	pred$dirDeathsPer1000personYr<-(pred$dirDeaths/(2*pred$nhost))*1000 
	return(pred)
}

# Function extractIndMort -extract INDIRECT MORTALITY estimates
extractIndMort<-function(scenNum,runDir, modelResultID) {
	# read in predictions for particular scenario number
	output<-read.table(paste(runDir,"wu_5000_", scenNum,".txt", sep="")) 
	# output is for 1 age groups 
	pred <- list()
	ageGroup <- 1
	# predicted result is all cause infant mortality rate
	pred$allCauseIMR <-output$V4[output$V2== ageGroup & output$V3==modelResultID$allCauseIMR]
	return(pred)
}

# =============================================================================
# OBECTIVES AND THEIR PLOT COMMANDS - functions to get field data and model predictions and to produce plots for each objective
# =============================================================================

#------------------------------------------------------------------------------
# Age Pattern of Incidence after intervention
#------------------------------------------------------------------------------
OBJ_AgePatternIncidenceAfterIntervention <- function(fieldDataAll, fieldDataID, modelResultID, densBias, runDir){
	# field data (scenario 30)
	incidenceField<-incidenceAfterInterventionFieldData(fieldDataAll, fieldDataID)
	# site names in incidenceField$site, scenario numbers in incidenceField$scenarioNum
	incidenceData <- list()
	incidenceData <- incidenceField$data
	
	# Incidence - model predictions	
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:1) {
		modelResults[[incidenceField$site[[scenario]]]] <- extractIncid(incidenceField$scenarioNum[scenario],runDir, modelResultID)
	}
		
	# LF calculation OBJ 1 - binomialLoglikelihood
	lossVectorObj_1 <- list()
	lossVectorObj_1$orig$total <- 0
	lossVectorObj_1$orig$valWithScenarioNum <- array(0,dim=c(length(incidenceField$scenarioNum),2))
  
	lossVectorObj_1$logLH$total <- 0
	lossVectorObj_1$logLH$valWithScenarioNum <- array(0,dim=c(length(incidenceField$scenarioNum),2))
	for (scenario in 1:1) {
		lossVectorObj_1$orig[[incidenceField$site[[scenario]]]] <-0
		lossVectorObj_1$logLH[[incidenceField$site[[scenario]]]] <-0
    
		for (agegroup in 1:length(incidenceData[[incidenceField$site[[scenario]]]]$age)){
			p <- modelResults[[incidenceField$site[[scenario]]]]$indivSurveys_prevalence[[agegroup]]
			k_mod <- modelResults[[incidenceField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			n_mod <- modelResults[[incidenceField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			k_mod_beta_one <- sum(k_mod)/sum(n_mod)*2
			k_mod_beta_two <- 2-k_mod_beta_one
			p[p>0.999]<-0.999
			p[p<0.001]<-0.001
			datanhost <- incidenceData[[incidenceField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			datanPositive  <- incidenceData[[incidenceField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			k_obs <- datanPositive
			n_obs <- datanhost 
			pData<-datanPositive/datanhost
			pData[pData>0.999]<-0.999
			pData[pData<0.001]<-0.001

			lossVectorObj_1$orig[[incidenceField$site[[scenario]]]] <- lossVectorObj_1$orig[[incidenceField$site[[scenario]]]]-sum(datanPositive*log ( p ) + ( datanhost-datanPositive ) *log( 1-p ))
			lossVectorObj_1$logLH[[incidenceField$site[[scenario]]]] <- lossVectorObj_1$logLH[[incidenceField$site[[scenario]]]] - sum( lbeta(k_obs+k_mod+k_mod_beta_one,(n_obs-k_obs)+(n_mod-k_mod)+k_mod_beta_two) - lbeta(k_mod+k_mod_beta_one,n_mod-k_mod+k_mod_beta_two) )
			rm(p, datanhost, datanPositive)
		}
		lossVectorObj_1$orig$valWithScenarioNum[scenario,] <- c(incidenceField$scenarioNum[scenario],lossVectorObj_1$orig[[incidenceField$site[[scenario]]]])
		lossVectorObj_1$logLH$valWithScenarioNum[scenario,] <- c(incidenceField$scenarioNum[scenario],lossVectorObj_1$logLH[[incidenceField$site[[scenario]]]])
		
		lossVectorObj_1$orig$total <- lossVectorObj_1$orig$total + lossVectorObj_1$orig[[incidenceField$site[[scenario]]]]
		lossVectorObj_1$logLH$total <- lossVectorObj_1$logLH$total + lossVectorObj_1$logLH[[incidenceField$site[[scenario]]]]
	}
	
	# PLOTS - INCIDENCE AFTER INTERVENTION
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(1,1,2,1))
	plot(incidenceData$Matsari$age,incidenceData$Matsari$prevalence,ylim=c(0,0.4),pch=15,log="x",xlim=c(0.1,100),ylab="prevalence",xlab="age")
	lines(incidenceData$Matsari$age, modelResults$Matsari$prevalence)
	legend("topright", c("field data","predictions"),lty=c(0,1),pch=c(15,NA),inset=.02,box.lty=0)
	mtext("Age patterns of prevalence of infection after intevention (Matsari)", NORTH<-3, line=0, adj=0.5, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_1)
}

#------------------------------------------------------------------------------
# Age Pattern of Prevalence And Parasite Density
#------------------------------------------------------------------------------
OBJ_AgePatternPrevalenceAndDensity<-function(fieldDataAll, fieldDataID, modelResultID, densBias, runDir){
	# field data (scenarios 24,28,29,35,34 and 31)
	asexualField<-asexualFieldData(fieldDataAll, fieldDataID)
	# site names in asexualField$site, scenario numbers in asexualField$scenarioNum
	asexualData <- list()
	asexualData <- asexualField$data
	
	# density bias vector, one for each scenario
	densBiasVec<-c(densBias$garki, densBias$garki, densBias$garki, densBias$nonGarki, densBias$nonGarki, densBias$nonGarki)
	# Asexual - model predictions	
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:6) {
		densBiasVal <- densBiasVec[scenario]
		modelResults[[asexualField$site[[scenario]]]] <- extractAS(asexualField$scenarioNum[scenario],runDir, densBiasVal, modelResultID)
	}
	
	# LF calculation OBJ 2 - binomialLoglikelihood
	lossVectorObj_2 <- list()
  
	lossVectorObj_2$orig$total <- 0
	lossVectorObj_2$orig$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))
  
	lossVectorObj_2$logLH$total <- 0
	lossVectorObj_2$logLH$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))
  
	for (scenario in 1:6) {
		lossVectorObj_2$orig[[asexualField$site[[scenario]]]] <-0
		lossVectorObj_2$logLH[[asexualField$site[[scenario]]]] <-0
		for (agegroup in 1:length(asexualData[[asexualField$site[[scenario]]]]$age)){
			p <- modelResults[[asexualField$site[[scenario]]]]$indivSurveys_prevalence[[agegroup]]
			k_mod <- modelResults[[asexualField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			n_mod <- modelResults[[asexualField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			k_mod_beta_one <- sum(k_mod)/sum(n_mod)*2
			k_mod_beta_two <- 2-k_mod_beta_one
			p[p>0.999]<-0.999
			p[p<0.001]<-0.001
			datanhost <- asexualData[[asexualField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			datanPositive  <- asexualData[[asexualField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			k_obs <- datanPositive
			n_obs <- datanhost
			lossVectorObj_2$orig[[asexualField$site[[scenario]]]] <- lossVectorObj_2$orig[[asexualField$site[[scenario]]]]-sum(datanPositive*log ( p ) + ( datanhost-datanPositive ) *log( 1-p ))
        pData<-datanPositive/datanhost
			  pData[pData>0.999]<-0.999
			  pData[pData<0.001]<-0.001
			  lossVectorObj_2$logLH[[asexualField$site[[scenario]]]] <- lossVectorObj_2$logLH[[asexualField$site[[scenario]]]] - sum( lbeta(k_obs+k_mod+k_mod_beta_one,(n_obs-k_obs)+(n_mod-k_mod)+k_mod_beta_two) - lbeta(k_mod+k_mod_beta_one,n_mod-k_mod+k_mod_beta_two) )
			  rm(p, datanhost, datanPositive)
      
		}
		
		lossVectorObj_2$orig$valWithScenarioNum[scenario,] <- c(asexualField$scenarioNum[scenario],lossVectorObj_2$orig[[asexualField$site[[scenario]]]])
		lossVectorObj_2$logLH$valWithScenarioNum[scenario,] <- c(asexualField$scenarioNum[scenario],lossVectorObj_2$logLH[[asexualField$site[[scenario]]]])
    
		lossVectorObj_2$orig$total <- lossVectorObj_2$orig$total + lossVectorObj_2$orig[[asexualField$site[[scenario]]]]
		lossVectorObj_2$logLH$total <- lossVectorObj_2$logLH$total + lossVectorObj_2$logLH[[asexualField$site[[scenario]]]]
    
	}
	
	# LF calculation OBJ 3 - log-normal log likelihood
	lossVectorObj_3 <- list()
  
	lossVectorObj_3$orig$total <- 0
	lossVectorObj_3$orig$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))
  
	lossVectorObj_3$logLH$total <- 0
	lossVectorObj_3$logLH$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))
  
	lossVectorObj_3$logLH_withoutConsts$total <- 0
	lossVectorObj_3$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))
  
	lossVectorObj_3$logLH_withConsts$total <- 0
	lossVectorObj_3$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(asexualField$scenarioNum),2))

	rho <- exp(-0.5*log(2*pi))
	for (scenario in 1:6) { 
		lossVectorObj_3$orig[[asexualField$site[[scenario]]]] <-0
		n<-0
		RSS <- 0

		for (agegroup in 1:length(asexualData[[asexualField$site[[scenario]]]]$age)){
			predV <- modelResults[[asexualField$site[[scenario]]]]$indivSurveys_sumlogDens[[agegroup]]/modelResults[[asexualField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			obsV <- asexualData[[asexualField$site[[scenario]]]]$indivSurveys_sumlogDens[[agegroup]]/asexualData[[asexualField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]
			RSS_notSummed <-(  predV-log(densBiasVec[scenario] )-obsV )^2
			RSS <- RSS + sum(RSS_notSummed[!is.na(RSS_notSummed)])
		
      #When Melissa introduced log likelihoods, the following was added to avoid log(0). Commenting out until we add log likelihoods back in (requires tackling TO DOs)  
      #if(any(modelResults[[asexualField$site[[scenario]]]]$indivSurveys_nPatent[[agegroup]]==0)){
			#	RSS <- NA
			#}
			n <- n + length(obsV[!is.na(obsV)])
			rm(predV, obsV)
		}

		# original loss function - residual sums of squares
		lossVectorObj_3$orig[[asexualField$site[[scenario]]]] <- RSS
		lossVectorObj_3$orig$valWithScenarioNum[scenario,] <- c(asexualField$scenarioNum[scenario], lossVectorObj_3$orig[[asexualField$site[[scenario]]]])
		lossVectorObj_3$orig$total <- lossVectorObj_3$orig$total + lossVectorObj_3$orig[[asexualField$site[[scenario]]]]

		sigma <-sqrt( RSS/(n-1)	)
		# terms of the log-likelihood can be removed as they are constant		
		lossVectorObj_3$logLH_withoutConsts[[asexualField$site[[scenario]]]] <- -(n*(-log(sigma+1)))
		# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
		lossVectorObj_3$logLH_withConsts[[asexualField$site[[scenario]]]] <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
		
		# we should really have log-likelihood with no constants
		# namely -(n*(-log(sigma))) BUT HAVE TO CHANGE WEIGHTS! and database
		# as 0.5*RSS/(sigma^2) = 0.5*(n-1)
		rm(RSS, n, sigma)

		lossVectorObj_3$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(asexualField$scenarioNum[scenario], lossVectorObj_3$logLH_withoutConsts[[asexualField$site[[scenario]]]])
		lossVectorObj_3$logLH_withoutConsts$total <- lossVectorObj_3$logLH_withoutConsts$total + lossVectorObj_3$logLH_withoutConsts[[asexualField$site[[scenario]]]]
    
		lossVectorObj_3$logLH_withConsts$valWithScenarioNum[scenario,] <- c(asexualField$scenarioNum[scenario], lossVectorObj_3$logLH_withConsts[[asexualField$site[[scenario]]]])
		lossVectorObj_3$logLH_withConsts$total <- lossVectorObj_3$logLH_withConsts$total + lossVectorObj_3$logLH_withConsts[[asexualField$site[[scenario]]]]

	}
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_3$orig <- lossVectorObj_3$logLH_withConsts
	lossVectorObj_3$logLH <- lossVectorObj_3$logLH_withoutConsts
	# LH - logliklihood (which is the same as orig loss Vector for this obj)
	# lossVectorObj_3$logLH$total <- lossVectorObj_3$orig$total_withConsts
	# lossVectorObj_3$logLH$valWithScenarioNum <- lossVectorObj_3$orig$valWithScenarioNum
	
		
	# PLOTS - PREVALENCE	
	par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), oma=c(2,2,2,2))
	textpos <- c(3,3,3,3,3,0.2)
	ylabels <- c("","","prevalence","","","")
	xlabels <- c("","","","","age","age")
	for (i in 1:6) {
		plot(asexualData[[i]]$age,asexualData[[i]]$prevalence,ylim=c(0,1),pch=15,ylab= ylabels[i],xlab= xlabels[i],log="x")
		lines(asexualData[[i]]$age, modelResults[[i]]$prevalence)	
		text(textpos[i],0.2,asexualData[[i]]$name)
	}
	mtext("Age patterns of prevalence of infection", NORTH<-3, line=0, adj=0.5, cex=1.2, col="black", outer=TRUE)
		
	# PLOTS - DENSITY	
	par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), oma=c(2,2,2,2))
	textpos <- c(3,3,3,3,3,0.2)
	ylabels <- c("","","density","","","")
	xlabels <- c("","","","","age","age")
	y_max <- c(10000,10000,10000,100000,100000,100000)
	for (i in 1:6) {
		plot(asexualData[[i]]$age,asexualData[[i]]$density,ylim=c(1,y_max[i]),pch=15,ylab= ylabels[i],xlab= xlabels[i],log="xy")
		lines(asexualData[[i]]$age,modelResults[[i]]$density)	
		text(textpos[i],5,asexualData[[i]]$name)
	}
	mtext("Age patterns of parasite density", NORTH<-3, line=0, adj=0.5, cex=1.2, col="black", outer=TRUE)
	
	combinedObj <- list()
	combinedObj$obj_2 <- lossVectorObj_2
	combinedObj$obj_3 <- lossVectorObj_3
	return(combinedObj)
}

#------------------------------------------------------------------------------
# Age Pattern of Number of concurrent infections
#------------------------------------------------------------------------------

OBJ_AgePatternMOI <- function(fieldDataAll, fieldDataID, modelResultID, runDir){	
	# field data (scenario 34)
	moiField<-moiFieldData(fieldDataAll, fieldDataID)
	# site names in moiField$site, scenario numbers in moiField$scenarioNum
	moiData <- list()
	moiData <- moiField$data
	
	# extract data for scenario 34 and calculate MOI (moiPred)
	# Asexual - model predictions	
	modelResults <- extractMoi(moiField$scenarioNum[1], runDir, modelResultID)

	# LF calculation OBJ 4 - poisson log-likelihood
	lossVectorObj_4 <- list()
	lossVectorObj_4$orig$total <- 0
	lossVectorObj_4$orig$valWithScenarioNum <- array(0,dim=c(length(moiField$scenarioNum),2))
	for (scenario in 1:1) { 
		lossVectorObj_4$orig[[moiField$site[[scenario]]]] <- 0
		lossVectorObj_4$orig$scen1_old_cpp <- 0
		for (agegroup in 1:length(moiData[[moiField$site[[scenario]]]]$age)){
			lambda <- modelResults$indivSurveys_moi[[agegroup]]*moiData[[moiField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			obs <- moiData[[moiField$site[[scenario]]]]$indivSurveys_numPatInf[[agegroup]]
      
      		# introduced by Ewan in 2017 see below, in 2020 April we removed due to testing
      		# alpha.est <- 1+modelResults$indivSurveys_numPatInf[[agegroup]]
      		# beta.est <- 1+modelResults$indivSurveys_nHosts[[agegroup]]
      		# p.est <- 1/(1+beta.est)
      		# r <- alpha.est*moiData[[moiField$site[[scenario]]]]$indivSurveys_nHosts[[agegroup]]
			
			tempVal <- -obs*log(obs/lambda)+obs-lambda
			tempVal[obs<=0.00001] <- -lambda[obs<=0.00001]
			lossVectorObj_4$orig[[moiField$site[[scenario]]]] <- lossVectorObj_4$orig[[moiField$site[[scenario]]]] - sum(tempVal)
			
			# introduced by Ewan in 2017 to take an approach to integrate out the uncertainty of the model predictions (from a finite number of model people) in each age/scenario bin with the corresponding conjugate statistical model.  For example, in the case of binomial data we can take the model positives, k_mod, and model population in that bin, n_mod, as  data to infer a posterior distribution for the model prevalence in that bin using the conjugate prior for the binomial which is the beta distribution, Beta(alpha,beta), giving a posterior, Beta(alpha+k_mod,beta+n_mod-k_mod).  [For choice of hyper-parameters alpha and beta I would cheat and use the mean prevalence across all age bins from the model, p_all_mod, weighted as two pseudo-observations: alpha = p_all_mod*2, beta = 2-alpha.]  The predictive distribution for the observed data in that bin is then the beta-binomial model giving a log-likelihood (ignoring constant terms) of log Beta(k_obs+k_mod+alpha,n_obs-k_obs+n_mod-k_mod+beta) - log Beta(k_mod+alpha,n_mod-k_mod+beta).  
			# lossVectorObj_4$orig[[moiField$site[[scenario]]]] <- lossVectorObj_4$orig[[moiField$site[[scenario]]]] - sum(dnbinom(obs,r,1-p.est,log=T))
			
			# this is the old calculation in cpp file, incorrect (TODO)
			tempVal_old_cpp <- -obs*log(obs/lambda)-obs+lambda			
			tempVal_old_cpp[obs<=0.00001] <- -lambda[obs<=0.00001]
			lossVectorObj_4$orig$scen1_old_cpp <- lossVectorObj_4$orig$scen1_old_cpp - sum(tempVal_old_cpp)
			rm(lambda, obs, tempVal, tempVal_old_cpp)
		}
		lossVectorObj_4$orig$valWithScenarioNum[scenario,] <- c(moiField$scenarioNum[scenario],lossVectorObj_4$orig[[moiField$site[[scenario]]]])
		lossVectorObj_4$orig$total <- lossVectorObj_4$orig$total + lossVectorObj_4$orig[[moiField$site[[scenario]]]]
		#lossVectorObj_4$orig$total <- lossVectorObj_4$orig$total + lossVectorObj_4$orig$scen1_old_cpp
	}
	# LH - logliklihood (which is the same as orig loss Vector for this obj)
	lossVectorObj_4$logLH$total <- lossVectorObj_4$orig$total
	lossVectorObj_4$logLH$valWithScenarioNum <- lossVectorObj_4$orig$valWithScenarioNum
		
	# PLOTS - MOI	
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(1,1,2,1))
	plot(moiData[[1]]$age,moiData[[1]]$moi,ylim=c(0,6),pch=15,log="x",xlim=c(0.1,100),ylab="MOI",xlab="age")
	lines(moiData[[1]]$age, modelResults$moi)
	legend("topright", c("field data","predictions"),lty=c(0,1),pch=c(15,NA),inset=.02,box.lty=0)
	mtext("Age pattern of number of concurrent infections", NORTH<-3, line=0, adj=0.5, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_4)
}

#------------------------------------------------------------------------------
# Age Pattern of Incidence of Clinical Malaria in Senegal
#------------------------------------------------------------------------------
OBJ_AgePatternIncidenceClinicalMalaria_Senegal <- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenarios 232 and 233)
	acuteSenegalField<-acuteSenegalFieldData(fieldDataAll, fieldDataID)
	# site names in acuteSenegalField$site, scenario numbers in acuteSenegalField$scenarioNum
	acuteSenegalData <- list()
	acuteSenegalData <- acuteSenegalField$data
	
	# acute - model predictions		
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:2) {
		# The rate multiplier is the duration in years for which episodes are collected
		# which is 5 years
		rateMultiplier<-5
		modelResults[[acuteSenegalField$site[[scenario]]]] <-extractAEpSenegal(acuteSenegalField$scenarioNum[scenario],runDir, rateMultiplier, modelResultID)
	}
	
	# TODO check survey 1 only	
	modelResults_survey1 <- list()
	for (scenario in 1:2) {
		# The rate multiplier is the duration in years for which episodes are collected
		# which is 5 years
		rateMultiplier<-5
		modelResults_survey1[[acuteSenegalField$site[[scenario]]]] <-extractAEpSenegal_survey1(acuteSenegalField$scenarioNum[scenario],runDir, rateMultiplier, modelResultID)
	}
	
	# LF calculation OBJ 5 - RSSbiased (non log) and log likelihood
	lossVectorObj_5a <- list()
	lossVectorObj_5a$orig$total <- 0
	lossVectorObj_5a$orig$valWithScenarioNum <- array(0,dim=c(length(acuteSenegalField$scenarioNum),2))
	lossVectorObj_5a$logLH$total <- 0
	lossVectorObj_5a$logLH$valWithScenarioNum <- array(0,dim=c(length(acuteSenegalField$scenarioNum),2))
	lossVectorObj_5a$logLH_withoutConsts$total <- 0
	lossVectorObj_5a$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(acuteSenegalField$scenarioNum),2))
	lossVectorObj_5a$logLH_withConsts$total <- 0
	lossVectorObj_5a$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(acuteSenegalField$scenarioNum),2))
	for (scenario in 1:2) {
		Residuals <- (acuteSenegalData[[acuteSenegalField$site[[scenario]]]]$incid- (modelResults[[acuteSenegalField$site[[scenario]]]]$incid))
		ResidSquared <- Residuals^2
		RSS <- sum(ResidSquared)
		# original loss function - residual sums of squares
		lossVectorObj_5a$orig[[acuteSenegalField$site[[scenario]]]] <- RSS
		# save values to loss vector - original
		lossVectorObj_5a$orig$valWithScenarioNum[scenario,] <- c(acuteSenegalField$scenarioNum[scenario], lossVectorObj_5a$orig[[acuteSenegalField$site[[scenario]]]])
		lossVectorObj_5a$orig$total <- lossVectorObj_5a$orig$total + lossVectorObj_5a$orig[[acuteSenegalField$site[[scenario]]]]
		# log-likelihood - loss function
		n <- length(ResidSquared)
		sigma <-sqrt( RSS/(n-1)	)		
		# terms of the log-likelihood can be removed as they are constant
		lossVectorObj_5a$logLH_withoutConsts[[acuteSenegalField$site[[scenario]]]] <- -(n*(-log(sigma+1)))
		# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
		lossVectorObj_5a$logLH_withConsts[[acuteSenegalField$site[[scenario]]]] <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
		rm(Residuals, ResidSquared, RSS, n, sigma)
		
		# save values to loss vector log-likelihood
		lossVectorObj_5a$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(acuteSenegalField$scenarioNum[scenario], lossVectorObj_5a$logLH_withoutConsts[[acuteSenegalField$site[[scenario]]]])
		lossVectorObj_5a$logLH_withoutConsts$total <- lossVectorObj_5a$logLH_withoutConsts$total + lossVectorObj_5a$logLH_withoutConsts[[acuteSenegalField$site[[scenario]]]]
		lossVectorObj_5a$logLH_withConsts$valWithScenarioNum[scenario,] <- c(acuteSenegalField$scenarioNum[scenario], lossVectorObj_5a$logLH_withConsts[[acuteSenegalField$site[[scenario]]]])
		lossVectorObj_5a$logLH_withConsts$total <- lossVectorObj_5a$logLH_withConsts$total + lossVectorObj_5a$logLH_withConsts[[acuteSenegalField$site[[scenario]]]]
	}
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_5a$logLH <- lossVectorObj_5a$logLH_withoutConsts
	
	# PLOT - age pattern of incidence of clinical malaria
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	plot(acuteSenegalData$Ndiop$age, modelResults$Ndiop$incid,ylim=c(0,7),log="x",type="l",xlab="age", ylab="episodes/py")
	lines(acuteSenegalData$Ndiop$age, acuteSenegalData$Ndiop$incid,type="b")
	lines(acuteSenegalData$Dielmo$age, modelResults$Dielmo$incid,lty=4)#type="l")
	lines(acuteSenegalData$Dielmo$age, acuteSenegalData$Dielmo$incid,type="b",pch=15)
	legend("topright", c("Ndiop field","Ndiop predictions","Dielmo field","Dielmo predictions"),lty=c(1,1,1,4),pch=c(21,NA,15,NA),inset=.02,box.lty=0)
	mtext("Age pattern of incidence of clinical malaria", NORTH<-3, line=1, adj=0.5, cex=1.2, col="black", outer=TRUE)
	mtext("(Ndiop and Dielmo:clinical malaria at health centre)", NORTH<-3, line=0, adj=0.65, cex=1.2, col="black", outer=TRUE)
	
	# # TEST TODO check survey 1 PLOT - age pattern of incidence of clinical malaria
	# par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	# plot(acuteSenegalData$Ndiop$age, modelResults_survey1$Ndiop$incid,ylim=c(0,7),log="x",type="l",xlab="age", ylab="episodes/py")
	# lines(acuteSenegalData$Ndiop$age, acuteSenegalData$Ndiop$incid,type="b")
	# lines(acuteSenegalData$Dielmo$age, modelResults_survey1$Dielmo$incid,lty=4)#type="l")
	# lines(acuteSenegalData$Dielmo$age, acuteSenegalData$Dielmo$incid,type="b",pch=15)
	# legend("topright", c("Ndiop field","Ndiop predictions","Dielmo field","Dielmo predictions"),lty=c(1,1,1,4),pch=c(21,NA,15,NA),inset=.02,box.lty=0)
	# mtext("TEST TEST SURVEY 1 to see difference ", NORTH<-3, line=1, adj=0.5, cex=1.9, col="black", outer=TRUE)
	# mtext("(Ndiop and Dielmo:clinical malaria at health centre)", NORTH<-3, line=0, adj=0.65, cex=1.2, col="black", outer=TRUE)

	return(lossVectorObj_5a)
}

#------------------------------------------------------------------------------
# Age Pattern of Incidence of Clinical Malaria in Idete 
#------------------------------------------------------------------------------

OBJ_AgePatternIncidenceClinicalMalaria_Idete <- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenario 49)
	acuteIdeteField <- acuteIdeteFieldData(fieldDataAll, fieldDataID)
	# site names in acuteIdeteField$site, scenario numbers in acuteIdeteField$scenarioNum
	acuteIdeteData <- list()
	acuteIdeteData <- acuteIdeteField$data
	
	# acute IDETE - model predictions		
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:1) {
		# The rate multiplier is 1/access, access = 36%
		accessConst<-2.7975236259999985
		rateMultiplier<-accessConst
		modelResults[[acuteIdeteField$site[[scenario]]]] <-extractAEpIdete(acuteIdeteField$scenarioNum[scenario],runDir, rateMultiplier, modelResultID)
	}

	# LF calculation OBJ 5 - RSSbiased (non log)
	lossVectorObj_5b <- list()
	lossVectorObj_5b$orig$total <- 0
	lossVectorObj_5b$orig$valWithScenarioNum <- array(0,dim=c(length(acuteIdeteField$scenarioNum),2))
	lossVectorObj_5b$logLH$total <- 0
	lossVectorObj_5b$logLH$valWithScenarioNum <- array(0,dim=c(length(acuteIdeteField $scenarioNum),2))
	lossVectorObj_5b$logLH_withoutConsts$total <- 0
	lossVectorObj_5b$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(acuteIdeteField $scenarioNum),2))
	lossVectorObj_5b$logLH_withConsts$total <- 0
	lossVectorObj_5b$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(acuteIdeteField $scenarioNum),2))
	for (scenario in 1:1) {
		# original loss function - residual sums of squares
		Residuals <- acuteIdeteData[[acuteIdeteField$site[[scenario]]]]$incid - (modelResults[[acuteIdeteField$site[[scenario]]]]$incid)
		ResidSquared <- Residuals^2
		RSS <- sum(ResidSquared)
		# save values to loss vector - original
		lossVectorObj_5b$orig[[acuteIdeteField$site[[scenario]]]] <- RSS
		lossVectorObj_5b$orig$valWithScenarioNum[scenario,] <- c(acuteIdeteField$scenarioNum[scenario], lossVectorObj_5b$orig[[acuteIdeteField$site[[scenario]]]])
		lossVectorObj_5b$orig$total <- lossVectorObj_5b$orig$total + lossVectorObj_5b$orig[[acuteIdeteField$site[[scenario]]]]
		# log-likelihood - loss function
		n <- length(ResidSquared)
		sigma <-sqrt( RSS/(n-1))
		# terms of the log-likelihood can be removed as they are constant# TODO HERE CHECK model 1001 gives negative :()
		lossVectorObj_5b$logLH_withoutConsts[[acuteIdeteField $site[[scenario]]]] <- -(n*(-log(sigma+1)))
		# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
		lossVectorObj_5b$logLH_withConsts[[acuteIdeteField $site[[scenario]]]] <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
		# TODO this is working for loglikelihood that is negative!!! TODO to fix
		# -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(RSS/sigma^2))
		#n*log(sigma)+0.5*(RSS/sigma^2)
		# to ensure we don't have a negative nlh we requre sigma > 1/[(2*pi)*exp(0.5*(1-n)/n)]
		# or RSS > (n-1)/[(2*pi)*exp((1-n)/n)]
		
		#nval<-40
		#(nval-1)/(2*pi)*exp((1-nval)/nval)# RSS
		#1/(2*pi)*exp(0.5*(1-nval)/nval) #sigma
		
		rm(Residuals, ResidSquared, RSS, n, sigma)
		# save values to loss vector log-likelihood
		lossVectorObj_5b$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(acuteIdeteField$scenarioNum[scenario], lossVectorObj_5b$logLH_withoutConsts[[acuteIdeteField$site[[scenario]]]])
		lossVectorObj_5b$logLH_withoutConsts$total <- lossVectorObj_5b$logLH_withoutConsts$total + lossVectorObj_5b$logLH_withoutConsts[[acuteIdeteField$site[[scenario]]]]
		lossVectorObj_5b$logLH_withConsts$valWithScenarioNum[scenario,] <- c(acuteIdeteField$scenarioNum[scenario], lossVectorObj_5b$logLH_withConsts[[acuteIdeteField$site[[scenario]]]])
		lossVectorObj_5b$logLH_withConsts$total <- lossVectorObj_5b$logLH_withConsts$total + lossVectorObj_5b$logLH_withConsts[[acuteIdeteField$site[[scenario]]]]
	}
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_5b$logLH <- lossVectorObj_5b$logLH_withoutConsts
	
	# PLOT - age pattern of incidence of clinical malaria - Idete
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	plot(acuteIdeteData$Idete$age, modelResults$Idete$incid,ylim=c(0,7),type="o",pch=22,xlab="age",ylab="episodes/py")
	lines(acuteIdeteData$Idete$age,modelResults$Idete$incid*accessConst,type="l",pch=0)
	lines(acuteIdeteData$Idete$age, acuteIdeteData$Idete$incid,type="b",pch=15)
	legend("topright", c("predicted total incidence","health centre data","36% predicted"), lty=c(1,1,1), pch=c(-1,15,0), merge=TRUE, inset=.02, box.lty=0)
	mtext("Age pattern of incidence of clinical malaria", NORTH<-3, line=1, adj=0.5, cex=1.2, col="black", outer=TRUE)
	mtext("(Idete predicted acute incidence)", NORTH<-3, line=0, adj=0.38, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_5b)
}

#------------------------------------------------------------------------------
# Age Pattern of parasite density threshold for clinical attack
#------------------------------------------------------------------------------
OBJ_AgePatternThresholdClinicalAttack <- function(fieldDataAll, fieldDataID, modelResultID, densBias, runDir){
	# field data (scenario 234)
	acutePyGenThresField <- acutePyGenThresFieldData(fieldDataAll, fieldDataID)
	# site names in acutePyGenThresField$site, scenario numbers in acutePyGenThresField$scenarioNum
	acutePyGenThresData <- list()
	acutePyGenThresData <- acutePyGenThresField$data
	# field data (from NdiopDielmo.csv)
	
	# acute IDETE - model predictions		
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:1) {
		# The rate multiplier associated with log parasite/leucocyt ratio and thus μ = 1/(8000densbias) is the non-garki density bias.
		biasPy <-1/(8000*densBias$nonGarki) # have checked this is correct 1/(8000*densBias$nonGarki)
		modelResults[[acutePyGenThresField$site[[scenario]]]] <-extractPyGen(acutePyGenThresField$scenarioNum[scenario],runDir, biasPy, modelResultID)
	}
	
	# LF calculation OBJ 6 - RSSbiased (log)
	lossVectorObj_6 <- list()
	lossVectorObj_6$orig$total <- 0
	lossVectorObj_6$orig$valWithScenarioNum <- array(0,dim=c(length(acutePyGenThresField$scenarioNum),2))
	lossVectorObj_6$logLH$total <- 0
	lossVectorObj_6$logLH$valWithScenarioNum <- array(0,dim=c(length(acutePyGenThresField$scenarioNum),2))
	lossVectorObj_6$logLH_withoutConsts$total <- 0
	lossVectorObj_6$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(acutePyGenThresField$scenarioNum),2))
	lossVectorObj_6$logLH_withConsts$total <- 0
	lossVectorObj_6$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(acutePyGenThresField$scenarioNum),2))
	for (scenario in 1:1) {
		#RSS_1 <- log(acutePyGenThresData[[acutePyGenThresField$site[[scenario]]]]$Ystar) - log(modelResults[[acutePyGenThresField$site[[scenario]]]]$pyrogt)
		# this is to check cpp if the use of biasPy is correct
		Residuals <- log(acutePyGenThresData[[acutePyGenThresField$site[[scenario]]]]$Ystar) - (modelResults[[acutePyGenThresField$site[[scenario]]]]$logpyrogt/modelResults[[acutePyGenThresField$site[[scenario]]]]$nhost) -log(biasPy)
		Residuals[is.na(acutePyGenThresData[[acutePyGenThresField$site[[scenario]]]]$Ystar)] <- 0
		ResidSquared <- Residuals^2
		RSS <- sum(ResidSquared)
		# original loss function - residual sums of squares
		lossVectorObj_6$orig[[acutePyGenThresField$site[[scenario]]]] <- RSS
		lossVectorObj_6$orig$valWithScenarioNum[scenario,] <- c(acutePyGenThresField$scenarioNum[scenario], lossVectorObj_6$orig[[acutePyGenThresField$site[[scenario]]]])
		lossVectorObj_6$orig$total <- lossVectorObj_6$orig$total + lossVectorObj_6$orig[[acutePyGenThresField$site[[scenario]]]]
		# log-likelihood - loss function
		n <- length(ResidSquared)
		sigma <-sqrt( RSS/(n-1)	)		
		# terms of the log-likelihood can be removed as they are constant
		lossVectorObj_6$logLH_withoutConsts[[acutePyGenThresField$site[[scenario]]]] <- -(n*(-log(sigma+1)))
		# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
		lossVectorObj_6$logLH_withConsts[[acutePyGenThresField$site[[scenario]]]] <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
		rm(Residuals, ResidSquared , RSS, n, sigma)
		# save values to loss vector log-likelihood
		lossVectorObj_6$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(acutePyGenThresField$scenarioNum[scenario], lossVectorObj_6$logLH_withoutConsts[[acutePyGenThresField$site[[scenario]]]])
		lossVectorObj_6$logLH_withoutConsts$total <- lossVectorObj_6$logLH_withoutConsts$total + lossVectorObj_6$logLH_withoutConsts[[acutePyGenThresField$site[[scenario]]]]
		lossVectorObj_6$logLH_withConsts$valWithScenarioNum[scenario,] <- c(acutePyGenThresField$scenarioNum[scenario], lossVectorObj_6$logLH_withConsts[[acutePyGenThresField$site[[scenario]]]])
		lossVectorObj_6$logLH_withConsts$total <- lossVectorObj_6$logLH_withConsts$total + lossVectorObj_6$logLH_withConsts[[acutePyGenThresField$site[[scenario]]]]	
	}
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_6$logLH <- lossVectorObj_6$logLH_withoutConsts
	
	# PLOT - age pattern of parasite density threshold for clinical attacks - Dielmo
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	plot(acutePyGenThresData$NdiopDielmo$age, acutePyGenThresData$NdiopDielmo$Ystar, pch=16, type="p", log="xy", ylim=c(0.1,5), xlim=c(0.1,100), ylab="parasite:leucocyte ratio", xlab="age")
	lines(acutePyGenThresData$NdiopDielmo$age, modelResults$NdiopDielmo$pyrogt)
	legend("topright", c("field data", "predictions"), lty=c(0,1), pch=c(16,NA), merge=TRUE, inset=.02, box.lty=0)
	mtext("Age pattern of parasite density threshold for clinical attacks", NORTH<-3, line=1, adj=0.5, cex=1.2, col="black", outer=TRUE)
	mtext("(Dielmo pyrogenic threshold)", NORTH<-3, line=0, adj=0.2, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_6)
}

#------------------------------------------------------------------------------
# Severe episodes vs prevalence
#------------------------------------------------------------------------------
OBJ_SevereEpisodesVsPrevalence <- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)
	severeData <- severeMarshSnowFieldData(fieldDataAll, fieldDataID)
	
	# severe - model predictions		
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:length(severeData$scenarioNum)) {
		pred_temp <- extractSM(severeData$scenarioNum[scenario],runDir, modelResultID)
		modelResults$prevalence[scenario] <- pred_temp$prevalence
		modelResults$sevEpPer1000personYr[scenario] <- pred_temp$sevEpPer1000personYr
		rm(pred_temp)
	}	
	
	# LF calculation OBJ 7 - squaredDeviationHospitalSevereMonthly
	lossVectorObj_7 <- list()
	lossVectorObj_7$orig$total <- 0
	lossVectorObj_7$orig$valWithScenarioNum <- array(0,dim=c(length(severeData$scenarioNum),2))
	lossVectorObj_7$logLH$total <- 0
	lossVectorObj_7$logLH$valWithScenarioNum <- array(0,dim=c(length(severeData$scenarioNum),2))
	lossVectorObj_7$logLH_withoutConsts$total <- 0
	lossVectorObj_7$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(severeData$scenarioNum),2))
	lossVectorObj_7$logLH_withConsts$total <- 0
	lossVectorObj_7$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(severeData$scenarioNum),2))
	# access for presentation of severe cases for this scenario (fit in the base model. Fixed unless fit for at later date (possible TODO))
	AccessForSevereCases <- 0.48
	ResidSquared <- array(0,dim = c(length(severeData$scenarioNum),1))
	for (scenario in 1:length(severeData$scenarioNum)) {
		indexT <- 1
    	while ( severeData$prevalence[indexT] <= modelResults$prevalence[scenario] ) {
        	indexT <- indexT +1
        }
        
		a0<- severeData$prevalence[indexT-1]
    a1<- severeData$prevalence[indexT]
    	# multiply here by access as match cpp (as data already divided by access)
    	f0<- severeData$MarshSnow_episodes[indexT-1]*AccessForSevereCases
    	f1<- severeData$MarshSnow_episodes[indexT]*AccessForSevereCases
    	interpRate <- (( modelResults$prevalence[scenario]-a0 ) / ( a1-a0 ) )* ( f1-f0 ) +f0
    	rm(indexT, a0, a1, f0, f1)
		ResidSquared[scenario] <- (log( (AccessForSevereCases*modelResults$sevEpPer1000personYr[scenario])/interpRate) )^2
		#(log( (AccessForSevereCases*modelResults$sevEpPer1000personYr[scenario]))-log(interpRate) )^2
		#(log( (AccessForSevereCases*modelResults$sevEpPer1000personYr[scenario])/interpRate) )^2
		# original loss function - per scenario			
		lossVectorObj_7$orig$valWithScenarioNum[scenario,] <- c(severeData$scenarioNum[scenario], ResidSquared[scenario])
		# log-likelihood - loss function - unable to save per scenario as only one data point per scenario
		lossVectorObj_7$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(severeData$scenarioNum[scenario], NA)
		lossVectorObj_7$logLH_withConsts$valWithScenarioNum[scenario,] <- c(severeData$scenarioNum[scenario], NA)
	}
	# original loss function - residual sums of squares
	RSS <- sum(ResidSquared)
	lossVectorObj_7$orig$total <- RSS
	# log-likelihood - loss function
	n <- length(ResidSquared)
	sigma <-sqrt(RSS/(n-1)	)		
	# terms of the log-likelihood can be removed as they are constant
	lossVectorObj_7$logLH_withoutConsts$total<- -(n*(-log(sigma+1)))
	# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
	lossVectorObj_7$logLH_withConsts$total <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
	
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_7$logLH <- lossVectorObj_7$logLH_withoutConsts	
		
	# PLOT - (hospitalisation rates) severe episodes per 100 person year versus prevalence
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	plot(severeData$prevalence, severeData$MarshSnow_episodes,ylim=c(0,100),type="o",pch=16,xlab="prevalence in children 0-9 years", ylab="episodes/1000 pyar in children 0-9y",)
	# add predictions
	points(modelResults$prevalence, modelResults$sevEpPer1000personYr)
	legend("topright", c("field data","predictions"),lty=c(1,0), pch=c(16,1),inset=.02,box.lty=0)
	mtext("Hospitalisation rate in relation to prevalence in children", NORTH<-3, line=1, adj=0.6, cex=1.2, col="black", outer=TRUE)
	mtext("(severe malaria episodes)", NORTH<-3, line=0, adj=0.25, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_7)
}

#------------------------------------------------------------------------------
# Age Pattern of Severe hospitalisation (0-9 years)
#------------------------------------------------------------------------------
OBJ_AgePatternOfSevere<- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenarios 158,167,173,176)
	severeField <- severeSnow1997FieldData(fieldDataAll, fieldDataID)
	# site names in severeField$site, scenario numbers in severeField$scenarioNum
	severeData <- list()
	severeData <- severeField$data
	
	# Severe - model predictions	
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:length(severeField$scenarioNum)) {
		modelResults[[severeField$site[[scenario]]]] <-extractSMRR(severeField$scenarioNum[scenario],runDir, modelResultID)
	}
	
	# LF calculation OBJ 8 - RSSHospitalSevereRR
	lossVectorObj_8 <- list()
	lossVectorObj_8$orig$total <- 0
	lossVectorObj_8$orig$valWithScenarioNum <- array(0,dim=c(length(severeField$scenarioNum),2))
	lossVectorObj_8$logLH$total <- 0
	lossVectorObj_8$logLH$valWithScenarioNum <- array(0,dim=c(length(severeField$scenarioNum),2))
	lossVectorObj_8$logLH_withConsts$total <- 0
	lossVectorObj_8$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(severeField$scenarioNum),2))
	lossVectorObj_8$logLH_withoutConsts$total <- 0
	lossVectorObj_8$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(severeField$scenarioNum),2))
  
	plotData <- list()
	plotData$orig <- array(0,dim=c(length(severeField$scenarioNum),2))
	plotData$pred <- array(0,dim=c(length(severeField$scenarioNum),2))
  
  
	for (scenario in 1:length(severeField$scenarioNum)) {
		predRR <- array(0,dim<-c(2,1))
		RR <- array(0,dim<-c(2,1))
		predRR[1] <- modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[2] / modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[1]
		predRR[2] <- modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[3] / modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[1]
		RR[1] <- (log( (predRR[1])/severeData[[severeField$site[[scenario]]]]$relativeRiskComparedAgeGroup1[2]) )
		RR[2] <- (log( (predRR[2])/severeData[[severeField$site[[scenario]]]]$relativeRiskComparedAgeGroup1[3]) )
    
		plotData$orig[[scenario,1]] <- severeData[[severeField$site[[scenario]]]]$relativeRiskComparedAgeGroup1[2]
		plotData$orig[[scenario,2]] <- severeData[[severeField$site[[scenario]]]]$relativeRiskComparedAgeGroup1[3]
		plotData$pred[[scenario,1]] <- predRR[1]
		plotData$pred[[scenario,2]] <- predRR[2]
    
		# original loss function - residual sums of squares
		lossVectorObj_8$orig[[severeField$site[[scenario]]]] <- sum(RR^2)
		lossVectorObj_8$orig$valWithScenarioNum[scenario,] <- c(severeField$scenarioNum[scenario], lossVectorObj_8$orig[[severeField$site[[scenario]]]])
		lossVectorObj_8$orig$total <- lossVectorObj_8$orig$total + lossVectorObj_8$orig[[severeField$site[[scenario]]]]
		rm(predRR,RR)
		# log-likelihood - loss function - unable to save per scenario as only one data point per scenario
		lossVectorObj_8$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(severeField$scenarioNum[scenario], NA)
		lossVectorObj_8$logLH_withConsts$valWithScenarioNum[scenario,] <- c(severeField$scenarioNum[scenario], NA)
	}
	# log-likelihood - loss function
	n <- length(lossVectorObj_8$orig$valWithScenarioNum[,2])
	RSS <- sum(lossVectorObj_8$orig$valWithScenarioNum[,2])
	sigma <-sqrt(RSS/(n-1)	)		
	# terms of the log-likelihood can be removed as they are constant
	# lossVectorObj_8$logLH_withoutConsts$total<- -(n*(-log(sigma)))
	# add the +1 in 2017
	lossVectorObj_8$logLH_withoutConsts$total<- -(n*(-log(sigma+1)))
	
	# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
	lossVectorObj_8$logLH_withConsts$total <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	# lossVectorObj_8$logLH <- lossVectorObj_8$logLH_withConsts
	# add 2017
	lossVectorObj_8$logLH <- lossVectorObj_8$logLH_withoutConsts
	
	# PLOT - MALARIA SEVERE RELATIVE RISK
	par(mfrow=c(3,1), mar=c(2,4,0.5,0.5), oma=c(2,2,2,2))
	# data
	scenario <-1
	plot(severeData[[severeField$site[[scenario]]]]$age,severeData[[severeField$site[[scenario]]]]$rate_epPer1000perYr,ylim=c(0,200),xlab = " ", ylab="episodes per 1000py",xaxt="n",type="b")
	for (scenario in 2:length(severeField$scenarioNum)) {
		lines(severeData[[severeField$site[[scenario]]]]$age,severeData[[severeField$site[[scenario]]]]$rate_epPer1000perYr,type="b")
	}
	axis(1,at=c(0.5,3,7),labels=c("1-11 months","1-4 years","5-9 years"))
	legend("topright", c("observed data (Snow et al 1997)"),lty=c(1,0), pch=c(21,1),inset=.02,box.lty=0)
	# model predictions
	scenario <-1
	plot(severeData[[severeField$site[[scenario]]]]$age, modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[1:3],ylim=c(0,200),pch=16,type="b",xlab = "age",ylab="episodes per 1000py",xaxt="n")	
	for (scenario in 2:length(severeField$scenarioNum)) {
		lines(severeData[[severeField$site[[scenario]]]]$age, modelResults[[severeField$site[[scenario]]]]$sevEpPer1000personYr[1:3],type="b",pch=16) # only first three age groups
	}
	axis(1,at=c(0.5,3,7),labels=c("1-11 months","1-4 years","5-9 years"))
	legend("topright", c("model predictions"),lty=c(1,0), pch=c(16,1),inset=.02,box.lty=0)
	# plot title
	mtext("Age pattern of hospitalisation : severe malaria", NORTH<-3, line=0, adj=0.8, cex=1.2, col="black", outer=TRUE)
  
	# PLOT - MALARIA SEVERE RELATIVE RISK SLOPE
	scenario <-1
	# means observed data
	plot(severeField$scenarioNum, plotData$orig[,1], ylim=c(0,2),xlim=c(0,4),pch=15,xlab="scenario", ylab="slope",)
	for (scenario in 1:length(severeField$scenarioNum)) {
	  points(scenario, plotData$orig[[scenario]],pch=15)
	}
	# model predictions
	for (scenario in 1:length(severeField$scenarioNum)) {
	  points(scenario, plotData$pred[[scenario]],pch=1)
	}
	legend("topright", c("field data","predictions"),pch=c(15,1),inset=.02,box.lty=0)
	
	
	return(lossVectorObj_8)
}

#------------------------------------------------------------------------------
# Direct Malaria Mortality
#------------------------------------------------------------------------------
OBJ_DirectMalariaMortality <- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenarios 301,302,303,312,316,317,318,326,327)
	dirMortField <- directMortalityFieldData(fieldDataAll, fieldDataID)
	# site names in dirMortField$site, scenario numbers in dirMortField$scenarioNum
	dirMortData <- list()
	dirMortData <- dirMortField$data

	# direct mortality - model predictions	
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:length(dirMortField$scenarioNum)) {
		modelResults[[dirMortField$site[[scenario]]]] <-extractMort(dirMortField$scenarioNum[scenario],runDir, modelResultID)
	}
	
	# LF calculation OBJ 9 - squaredDeviationlogRateMM
	lossVectorObj_9 <- list()
	lossVectorObj_9$orig$total <- 0
	lossVectorObj_9$orig$valWithScenarioNum <- array(0,dim=c(length(dirMortField$scenarioNum),2))
	lossVectorObj_9$logLH$total <- 0
	lossVectorObj_9$logLH$valWithScenarioNum <- array(0,dim=c(length(dirMortField$scenarioNum),2))
	lossVectorObj_9$logLH_withConsts$total <- 0
	lossVectorObj_9$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(dirMortField$scenarioNum),2))
	lossVectorObj_9$logLH_withoutConsts$total <- 0
	lossVectorObj_9$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(dirMortField$scenarioNum),2))
	for (scenario in 1:length(dirMortField$scenarioNum)) {
		pred <- modelResults[[dirMortField$site[[scenario]]]]$dirDeaths/(2*modelResults[[dirMortField$site[[scenario]]]]$nhost)
		Residuals <- log( pred/dirMortData[[dirMortField$site[[scenario]]]]$dirDeathRateForAgeGroup1)
		ResidSquared <- Residuals^2
		RSS <- sum(ResidSquared)
		lossVectorObj_9$orig[[dirMortField$site[[scenario]]]] <- RSS
		lossVectorObj_9$orig$valWithScenarioNum[scenario,] <- c(dirMortField$scenarioNum[scenario], lossVectorObj_9$orig[[dirMortField$site[[scenario]]]])
		lossVectorObj_9$orig$total <- lossVectorObj_9$orig$total + lossVectorObj_9$orig[[dirMortField$site[[scenario]]]]
		rm(Residuals, ResidSquared,RSS,pred)
		# log-likelihood - loss function - unable to save per scenario as only one data point per scenario
		lossVectorObj_9$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(dirMortField$scenarioNum[scenario], NA)
		lossVectorObj_9$logLH_withConsts$valWithScenarioNum[scenario,] <- c(dirMortField$scenarioNum[scenario], NA)
	}
	# log-likelihood - loss function
	n <- length(lossVectorObj_9$orig$valWithScenarioNum[,2])
	RSS <- sum(lossVectorObj_9$orig$valWithScenarioNum[,2])
	sigma <-sqrt(RSS/(n-1)	)		
	# terms of the log-likelihood can be removed as they are constant
	lossVectorObj_9$logLH_withoutConsts$total<- -(n*(-log(sigma+1)))
	# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
	lossVectorObj_9$logLH_withConsts$total <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_9$logLH <- lossVectorObj_9$logLH_withoutConsts	
	
	# PLOT - under 5 direct malaria mortality
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	# data
	scenario <-1
	# means observed data
	plot(dirMortData[[dirMortField$site[[scenario]]]]$EIR, dirMortData[[dirMortField$site[[scenario]]]]$dirDeathRate_epPer1000perYr, ylim=c(0,40),xlim=c(5,400),log="x",pch=15,xlab="EIR", ylab="deaths/1000 person year",)
	for (scenario in 2:length(dirMortField$scenarioNum)) {
		points(dirMortData[[dirMortField$site[[scenario]]]]$EIR, dirMortData[[dirMortField$site[[scenario]]]]$dirDeathRate_epPer1000perYr,pch=15)
	}
	# max and min observed data
	for (scenario in 1:length(dirMortField$scenarioNum)) {
		arrows(dirMortData[[dirMortField$site[[scenario]]]]$EIR, dirMortData[[dirMortField$site[[scenario]]]]$dirDeathRate_epPer1000perYr_UpperCI, dirMortData[[dirMortField$site[[scenario]]]]$EIR, dirMortData[[dirMortField$site[[scenario]]]]$dirDeathRate_epPer1000perYr_LowerCI,code=3,length=.04,angle=90) # confidence intervals
	}
	# model predictions
	for (scenario in 1:length(dirMortField$scenarioNum)) {
		points(dirMortData[[dirMortField$site[[scenario]]]]$EIR, modelResults[[dirMortField$site[[scenario]]]]$dirDeathsPer1000personYr,pch=1)
	}
	 
	legend("topleft", c("field data","predictions"),pch=c(15,1),inset=.02,box.lty=0)
	mtext("Malaria specific mortality in children (Under 5y)", NORTH<-3, line=0, adj=0.6, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_9)
}

#------------------------------------------------------------------------------
# All cause mortality
#------------------------------------------------------------------------------
OBJ_indirectMortality <- function(fieldDataAll, fieldDataID, modelResultID, runDir){
	# field data (scenarios 401,402,408,411,414,415,416,417,418,422,426)
	indirMortField <- indirectMortalityFieldData(fieldDataAll, fieldDataID)
	# site names in indirMortField$site, scenario numbers in indirMortField$scenarioNum
	indirMortData <- list()
	indirMortData <- indirMortField$data
	
	# direct mortality - model predictions	
	# extract data for each of the scenarios and store to modelResults
	modelResults <- list()
	for (scenario in 1:length(indirMortField$scenarioNum)) {
		modelResults[[indirMortField$site[[scenario]]]] <-extractIndMort(indirMortField$scenarioNum[scenario],runDir, modelResultID)
	}
	
	# LF calculation OBJ 10 - squaredDeviationlogRateIM
	lossVectorObj_10 <- list()
	lossVectorObj_10$orig$total <- 0
	lossVectorObj_10$orig$valWithScenarioNum <- array(0,dim=c(length(indirMortField$scenarioNum),2))
	lossVectorObj_10$logLH$total <- 0
	lossVectorObj_10$logLH$valWithScenarioNum <- array(0,dim=c(length(indirMortField$scenarioNum),2))
	lossVectorObj_10$logLH_withConsts$total <- 0
	lossVectorObj_10$logLH_withConsts$valWithScenarioNum <- array(0,dim=c(length(indirMortField$scenarioNum),2))
	lossVectorObj_10$logLH_withoutConsts$total <- 0
	lossVectorObj_10$logLH_withoutConsts$valWithScenarioNum <- array(0,dim=c(length(indirMortField$scenarioNum),2))
	for (scenario in 1:length(indirMortField$scenarioNum)) {
		Residuals <- log( modelResults[[indirMortField$site[[scenario]]]]$allCauseIMR
/indirMortData[[indirMortField$site[[scenario]]]]$indirDeath)
		ResidSquared <- (Residuals)^2
		RSS <- sum(ResidSquared)
		lossVectorObj_10$orig[[indirMortField$site[[scenario]]]] <- RSS
		lossVectorObj_10$orig$valWithScenarioNum[scenario,] <- c(indirMortField$scenarioNum[scenario], lossVectorObj_10$orig[[indirMortField$site[[scenario]]]])
		lossVectorObj_10$orig$total <- lossVectorObj_10$orig$total + lossVectorObj_10$orig[[indirMortField$site[[scenario]]]]
		rm(Residuals, ResidSquared, RSS)
		# log-likelihood - loss function - unable to save per scenario as only one data point per scenario
		lossVectorObj_10$logLH_withoutConsts$valWithScenarioNum[scenario,] <- c(indirMortField$scenarioNum[scenario], NA)
		lossVectorObj_10$logLH_withConsts$valWithScenarioNum[scenario,] <- c(indirMortField$scenarioNum[scenario], NA)
	}
	# log-likelihood - loss function
	n <- length(lossVectorObj_10$orig$valWithScenarioNum[,2])
	RSS <- sum(lossVectorObj_10$orig$valWithScenarioNum[,2])
	sigma <-sqrt(RSS/(n-1)	)		
	# terms of the log-likelihood can be removed as they are constant
	lossVectorObj_10$logLH_withoutConsts$total<- -(n*(-log(sigma+1)))
	# note last term in equation is 0.5*sum(RSS)/(sigma^2)), which simplifies to 0.5*(n-1) given the definition of sigma
	lossVectorObj_10$logLH_withConsts$total <- -(n*((-0.5*log(2*pi))-log(sigma))-0.5*(n-1))
	# set the log-likelihood to be logLH with constants TODO decide what to do here, some models, logLH without constants is negative
	lossVectorObj_10$logLH <- lossVectorObj_10$logLH_withoutConsts	
	
	
	# PLOT - all causes infant indirect mortality
	par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=c(2,1,3,1))
	# data
	scenario <-1
	# observed data
	plot(indirMortData[[indirMortField$site[[scenario]]]]$EIR,  indirMortData[[indirMortField$site[[scenario]]]]$indirDeath,ylim=c(60,140),xlim=c(5,400),log="x",pch=15,xlab="EIR", ylab="deaths/1000 livebirths")
	for (scenario in 2:length(indirMortField$scenarioNum)) {
		points(indirMortData[[indirMortField$site[[scenario]]]]$EIR, indirMortData[[indirMortField$site[[scenario]]]]$indirDeath,pch=15)
	}
	# model predictions
	for (scenario in 1:length(indirMortField$scenarioNum)) {
		points(indirMortData[[indirMortField$site[[scenario]]]]$EIR, modelResults[[indirMortField$site[[scenario]]]]$allCauseIMR,pch=1)
	}
	legend("topleft", c("field data","predictions"),pch=c(15,1),inset=.02,box.lty=0)
	mtext("All-cause infant mortality rate", NORTH<-3, line=0, adj=0.6, cex=1.2, col="black", outer=TRUE)
	
	return(lossVectorObj_10)
}
	

# =============================================================================
# MAIN FUNCTION TO PLOT ALL RESULTS - allPlots
# =============================================================================
allPlots<-function(runID,dirname,params){
	runDir <- paste(sep = "",dirname,runID,"/")

        # name the pdf for plotting output
	pdf(file=paste(dirname,runID,".pdf", sep=""),paper="a4")
	
	# -----------------------------------------------------------------
	# read in all field data
	fieldDataAll <-read.table(paste("./fieldData",".txt", sep=""))
	
	# get codes for field data ids
	fieldDataID <- fieldData_nameCodes()
	
	# get codes for model Results ids
	modelResultID <- modelResults_nameCodes()
	# --------------------------------------------------------------------
	# density bias for observed parasite densities in garki sites and non-garki
        densBias <- list()
        densBias$nonGarki <- params$DensityBiasnonGarki
        densBias$garki <- params$DensityBiasGarki
	# --------------------------------------------------------------------
	# set loss Vector etc
	lossVector <- array(NA,dim=c(11,1))
	lossVector_logLH <- array(NA,dim=c(11,1))
	LF <- list()
	# --------------------------------------------------------------------
	# INCIDENCE after intervention - 
	# (scenarios 30)
	LF$obj_1 <- OBJ_AgePatternIncidenceAfterIntervention(fieldDataAll, fieldDataID, modelResultID, densBias, runDir)
	lossVector[1] <- LF$obj_1$orig$total
	lossVector_logLH[1] <- LF$obj_1$logLH$total
	# --------------------------------------------------------------------
	# ASEXUAL - prevalence and density (2 plots produced)
	# (scenarios 24,28,29,35,34 and 31)
	LF_temp <- list()
	LF_temp <- OBJ_AgePatternPrevalenceAndDensity(fieldDataAll, fieldDataID, modelResultID, densBias, runDir)
	LF$obj_2 <- LF_temp$obj_2
	LF$obj_3 <- LF_temp$obj_3
	rm(LF_temp)
	lossVector[2] <- LF$obj_2$orig$total
	lossVector[3] <- LF$obj_3$orig$total
	lossVector_logLH[2] <- LF$obj_2$logLH$total
	lossVector_logLH[3] <- LF$obj_3$logLH$total
	# --------------------------------------------------------------------
	# MOI - multiplicity of infection
	# (scenario 34)
	LF$obj_4 <- OBJ_AgePatternMOI(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[4] <- LF$obj_4$orig$total
	lossVector_logLH[4] <- LF$obj_4$logLH$total
	# --------------------------------------------------------------------
	# ACUTE EPISODES: NDIOP AND DIELMO INCIDENCE
	# (scenarios 232 and 233)
	LF$obj_5a <- OBJ_AgePatternIncidenceClinicalMalaria_Senegal(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[5] <- LF$obj_5a$orig$total
	lossVector_logLH[5] <- LF$obj_5a$logLH$total
	# --------------------------------------------------------------------
	# ACUTE EPISODES: IDETE
	# (scenario 49)
	LF$obj_5b <- OBJ_AgePatternIncidenceClinicalMalaria_Idete(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[6] <- LF$obj_5b$orig$total
	lossVector_logLH[6] <- LF$obj_5b$logLH$total
	# --------------------------------------------------------------------
	# ACUTE EPISODES: DIELMO PYROGENIC THRESHOLD
	# (scenario 234)
	LF$obj_6 <- OBJ_AgePatternThresholdClinicalAttack(fieldDataAll, fieldDataID, modelResultID, densBias, runDir)
	lossVector[7] <- LF$obj_6$orig$total
	lossVector_logLH[7] <- LF$obj_6$logLH$total
	# --------------------------------------------------------------------
	# SEVERE MALARIA (prevalence vs episodes)
	# (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)
	LF$obj_7 <- OBJ_SevereEpisodesVsPrevalence(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[8] <- LF$obj_7$orig$total
	lossVector_logLH[8] <- LF$obj_7$logLH$total
	# ---------------------------------------------------------------------
	# SEVERE MALARIA RR (relative risk)
	# (scenarios 158,167,173,176)
	LF$obj_8 <- OBJ_AgePatternOfSevere(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[9] <- LF$obj_8$orig$total
	lossVector_logLH[9] <- LF$obj_8$logLH$total
	# ---------------------------------------------------------------------
	# DIRECT MORTALITY
	# (scenarios 301,302,303,312,316,317,318,326,327)
	LF$obj_9 <- OBJ_DirectMalariaMortality(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[10] <- LF$obj_9$orig$total	
	lossVector_logLH[10] <- LF$obj_9$logLH$total
	# ---------------------------------------------------------------------
	# INDIRECT MORTALITY
	# (scenarios 401,402,408,411,414,415,416,417,418,422,426)
	LF$obj_10 <- OBJ_indirectMortality(fieldDataAll, fieldDataID, modelResultID, runDir)
	lossVector[11] <- LF$obj_10$orig$total
	lossVector_logLH[11] <- LF$obj_10$logLH$total
	# ---------------------------------------------------------------------
	# LossFunction weights # TODO this is hard coded at the moment
	lossFunc_weights <- array(c(0.001, 0.001, 0.01, 0.01, 1, 1, 1, 2, 2, 1, 10),dim=c(11,1))
	# ---------------------------------------------------------------------
	# WEIGHTED LOSS FUNCTION
	weighted_LF <- list()
	weighted_LF$vector <- lossFunc_weights*lossVector
	weighted_LF$vector_new <- lossFunc_weights*lossVector_logLH
	weighted_LF$total <- sum(weighted_LF$vector)
	weighted_LF$total_avgOverScen <- weighted_LF$total/61
  
	## Printing of lf into table in pdf 
	tableToPrint <- array(0,dim=c(11,4))
	tableToPrint[,1] <- lossVector
	tableToPrint[,2] <- lossFunc_weights
	tableToPrint[,3] <- weighted_LF$vector
	tableToPrint[,4] <- weighted_LF$vector_new
	
	# not the most elegant way to print to pdf but is quick fix
	plot.new()
	limUpper_y <- 50
	plot.window(xlim=c(0,1),ylim=c(0,limUpper_y)) 
	title(paste("Summary of loss function - Model",runID,sep=" "),cex.main=1.3) 
	
	rowNames_descObj<-c("Age prevalence  \nafter intervention", "Age pattern \nof prevalence", "Age pattern \nof parasite density", "Multiplicity of \nInfection", "Age pattern of clinical \nincidence : Senegal", "Age pattern of clinical \nincidence : Idete", "Clinical Threshold", "Severe episodes \nwith prevalence", "Age pattern of \nsevere episodes", "Direct mortality \nwith EIR", "All-cause mortality \nwith EIR")
	rowNames_obj<-c("obj_1", "obj_2", "obj_3", "obj_4", "obj_5a", "obj_5b", "obj_6", "obj_7", "obj_8", "obj_9", "obj_10")
	rowNames_losVec<-c("lossVec_1", "lossVec_2", "lossVec_3", "lossVec_4", "lossVec_5", "lossVec_6", " lossVec_7", "lossVec_8", "lossVec_9", "lossVec_10", "lossVec_11")
	
	text(0.0,(limUpper_y-1), adj=c(0,0),lab=paste("total weighted loss function     ", as.character(format(weighted_LF$total_avgOverScen, digits = 4, nsmall = 5)),"(averaged)",sep="  "))
	text(0.0,(limUpper_y-3), adj=c(0,0),lab=paste("total weighted loss function  ", as.character(format(weighted_LF$total, digits = 4, nsmall = 5)),sep=" "))
  
	plot.new()
	tableGrid <- data.frame(rowNames_descObj,round(lossVector, digits=6),as.character(lossFunc_weights),round(weighted_LF$vector, digits=4),round(weighted_LF$vector_new, digits=4))
	colnames(tableGrid)=c( "objective", "unweighted LF","weights","weighted LF","new weighted LF")
	grid.table(tableGrid)
  
	dev.off()
  	write.csv(tableGrid,paste0(dirname,runID,".csv"))
  
    resultsToReturn <- list()
	resultsToReturn$lossVector <- lossVector
	resultsToReturn$lossVector_logLH <- lossVector_logLH
	resultsToReturn$weighted_LF <- weighted_LF
	return(resultsToReturn)
}
