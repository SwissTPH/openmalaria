/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include "transmissionModel.h"
#include "simulation.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include "noVectorControl.h"
#include "vectorControl.h"
#include <math.h> 
#include "global.h" 
#include <gsl/gsl_vector.h> 
#include "infection.h" 
#include "intervention.h" 
#include <iostream>

TransmissionModel::TransmissionModel(){
  strcpy(fnametestentopar,"output_ento_para.txt");
}

//static (class) variables
const double TransmissionModel::susceptibility= 0.702;
const double TransmissionModel::totalInfectionrateVariance= 1.0;
const double TransmissionModel::min_EIR_mult= 0.01; 
const double TransmissionModel::agemin[nwtgrps] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60 };
const double TransmissionModel::agemax[nwtgrps] = { 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 19.99, 24.99, 29.99, 39.99, 49.99, 59.99, 60.99 };
const double TransmissionModel::wtprop[nwtgrps] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
const double TransmissionModel::bsa_prop[nwtgrps] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

extern double BaselineAvailabilityShapeParam;

void TransmissionModel::initEntoParameters() {
  EIRRotateAngle = M_PI/2;
  nspore=get_nspore();
  gamma_p=get_parameter(5);
  Sinf=1-exp(-get_parameter(1));
  Simm=get_parameter(3);
  Estar=get_parameter(2);
  Xstar_p=get_parameter(4);
  BaselineAvailabilityShapeParam=get_parameter(16);
  maxIntervals=maxDurIntPhaseEIR/interval;
  no = (int *) malloc(((intervalsPerYear))*sizeof(int));
  kappa = (double *) malloc(((intervalsPerYear))*sizeof(double));
  initialKappa = (double *) malloc(((intervalsPerYear))*sizeof(double));
  EIR = new double[intervalsPerYear];
  origEIR = (double *) malloc(((intervalsPerYear))*sizeof(double));
  ino = (int *) malloc(((maxIntervals))*sizeof(int));
  inoX = maxIntervals;
  intEIR = (double *) malloc(((maxIntervals))*sizeof(double));
  intEIRX = maxIntervals;
  
  //! Expected number of inoculations
  /// Product of measured EIR, susceptibility and length of time interval
  double gsi = 1.0;
  
  // FCEIR[] is the array of parameters of the Fourier approximation to the annual EIR
  // We set these here for now - with no if statement.
  // We will need to deal with this cleanly later.
  // We use the order, a0, a1, b1, a2, b2, ...
  // TODO: Move to XML.
  FCEIRX = 5;
  FCEIR = (double *) malloc((FCEIRX)*sizeof(double));
  FCEIR[0] = -0.926517;
  FCEIR[1] = -0.692164;
  FCEIR[2] = 0.002098;
  FCEIR[3] = 0.401189;
  FCEIR[4] = -0.375356;

  //! constant defining the constraint for the Gamma shape parameters
  /// Used for the case where availability is assumed gamma distributed
  double r_square_Gamma;
  /*
  r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
  r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
  */
  r_square_Gamma=0.649;
  //such that r_square_LogNormal =0.5
  
  //! constant defining the constraint for the log Normal variance
  /// Used for the case where availability is assumed log Normally distributed
  double r_square_LogNormal = log(1.0+r_square_Gamma);
  
  //TODO: Sanity check for sqrt and division by zero
  if ( isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
    InfectionrateShapeParam = (BaselineAvailabilityShapeParam+1.0) / (r_square_Gamma*BaselineAvailabilityShapeParam - 1.0);
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  else if(isOptionIncluded(modelVersion, lognormalMassAction) ||
          isOptionIncluded(modelVersion, lognormalMassActionPlusPreImm)) {
    InfectionrateShapeParam = sqrt(r_square_LogNormal - 1.86*pow(BaselineAvailabilityShapeParam, 2));
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  inputEIR();
  initAgeExposureConversion();
}



double TransmissionModel::getRelativeAvailability(double ageyrs) {
  //60 yrs is the last cutpoint in the human growth curve
  double limit = std::min (ageyrs, 60.0);
  int i=0;
  while(agemax[i] < limit)
    ++i;
  
  // the ageSpecificRelativeAvailability vector contains proportions of the
  // total bites received by a host of this age when competing with an adult 
  return ageSpecificRelativeAvailability[i];
}

double TransmissionModel::getWeight(double age){
  // NOTE: did effectively run for i = -1; i <= nages; ++i
  for (int i=0; i < nages; i++) {
    if (agemax[i] > age)
      return 120.0 * wtprop[i];
  }
  return 1.0;
}

void TransmissionModel::initAgeExposureConversion(){
  for (int i=0; i<nages; i++) {
    ageSpecificRelativeAvailability[i] = bsa_prop[i] / (1-bsa_prop[i]);
  }

  // The following code calculates average bites for children <6 years old
  // relative to adults. This is for analysis of Saradidi (Beier et al) data.

  // avbites_6 is the exposure of children < 6 years old relative to adults
  double avbites_6=0.0;
  
  for (int i = 0; agemin[i] <  6; ++i) {
    //agemin0 is the lower bound of the youngest age group used in this calculation
    double agemin0=agemin[i];
    if (agemax[i] > 0.5) {
        if (agemin0 < 0.5)
          agemin0=0.5;
    
      avbites_6 += ageSpecificRelativeAvailability[i] * (agemax[i]-agemin0)/5.5;
    }
  }
  biteratio_6 = avbites_6 / (1.0-avbites_6);
}


void TransmissionModel::inputEIR () {
  // Reads in the estimates of the EIR for each village and each day
  // and converts this into EIR estimates per five day period
  // assuming that the annual cycle repeated during the pre-intervention period

  // TODO: Properly include new simulationMode for entomological model. 
  // Depending on flags, this should either take in given EIR and smooth through
  // a DFT or directly take in Fourier coefficients and create an EIR over time.

  //initialise all the EIR arrays to 0
  if ( simulationMode !=  transientEIRknown) {
    for (int j=0;j<intervalsPerYear; j++) {
      EIR[j]=0.0;
      no[j]=0;
    }
  } else {
    for (int j=0;j<maxIntervals; j++) {
      intEIR[j]=0.0;
      ino[j]=0;
    }
  }
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR();
  double EIRdaily=get_eir_daily(0);
  for (int mpcday=0; EIRdaily != missing_value;) {
    updateEIR(mpcday, std::max(EIRdaily, minEIR));
    ++mpcday;
    EIRdaily = get_eir_daily(mpcday);
  }
  
  // Calculate total annual EIR 
  if ( simulationMode !=  transientEIRknown) {
    annualEIR=0.0;
    for (int j=0;j<intervalsPerYear; j++) {
      annualEIR += interval*EIR[j];
    }
  } else {
    annualEIR=-9.99;
  }
  
  // We copy EIR into origEIR
  for(int j=0; j<intervalsPerYear; j++){
    origEIR[j] = EIR[j];
  }

  // For now we assume that we can manipulate EIR depending on the value of FTSmoothEIR.
  // We are assuming that the simulation mode is not set to transientEIRknown.
  // This needs to be rewritten properly once we have introduced a new simulationMode.
  // If FTSmoothEIR is 0, we do nothing.
  // If FTSmoothEIR is 1, we smooth the EIR using the first 3 modes of the discrete
  //		Fourier Transform
  if(FTSmoothEIR==1){
# ifdef TransmissionModel_PrintOrigEIR
    PrintArray(fnametestentopar, "originalEIR", origEIR, intervalsPerYear);
# endif
  
    logDFTThreeModeSmooth(EIR, origEIR, intervalsPerYear, intervalsPerYear);
  }
  if(ifUseFC==1){
    calcInverseDFTExp(EIR, intervalsPerYear, FCEIR, FCEIRX);
# ifdef TransmissionModel_PrintEIRaIDFT
    PrintArray(fnametestentopar, "EIRafterIDFT", EIR, intervalsPerYear);
# endif
  }
  
  if(ifrotateEIR){
    rotateArray(EIR, intervalsPerYear, EIRRotateAngle);
  }
}

void TransmissionModel::updateEIR (int day, double EIRdaily) {
  // Processes each daily EIR estimate, allocating each day in turn to the
  // appropriate time period. EIRdaily is the value of the daily EIR read in
  // from the .XML file :

  // istep is the time period to which the day is assigned.  The result of the
  // division is automatically rounded down to the next integer
  int istep= 1 + (day-1) / interval;
  if ( simulationMode !=  transientEIRknown) {
    int i1 = modIntervalsPerYear(istep) - 1;
    no[i1]++;
    //EIR() is the arithmetic mean of the EIRs assigned to the 73 different recurring time points
    EIR[i1] = ((EIR[i1] * (no[i1]-1)) + EIRdaily) / no[i1];
  } else {
    int i1=istep - 1;
    ino[i1]++;
    intEIR[i1]= ((intEIR[i1] * (ino[i1]-1)) + EIRdaily) / ino[i1];
  }
}


void TransmissionModel::clearTransmissionModelParameters () {

  free(no);
  free(kappa);
  free(initialKappa);
  delete [] EIR;
  free(ino);
  free(intEIR);
  free(origEIR);
  free(FCEIR);

}

double TransmissionModel::averageEIR () {
  // Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
  double valaverageEIR=0.0;
  int i = 0;
  for (; get_eir_daily(i) !=  missing_value; ++i) {
    valaverageEIR += get_eir_daily(i);
  }
  valaverageEIR /= i;
  return valaverageEIR;
}

double TransmissionModel::getExpectedNumberOfInfections (Human& human, double age_adj_EIR) {
  double efficacy = human.getPEVEfficacy();
  double baseAvailToMos = human.getBaselineAvailabilityToMosquitoes();
  double normp;
  double survivalOfInoculum;
  double Infectionrate;
  //The age-adjusted EIR, possibly adjusted for bed nets.
  double effectiveEIR;
  //TODO: we have never used this ITN model.  It can be removed.
  if ( ITN) {
    effectiveEIR=age_adj_EIR*sqrt(Pu1/Pu0);
  } else {
    effectiveEIR=age_adj_EIR;
  }
  double ExpectedInfectionRate = effectiveEIR * baseAvailToMos * susceptibility * interval;
  double expectedNumInfections;
  if ( isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
    Infectionrate=(W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
    expectedNumInfections=Infectionrate;
  } else if( isOptionIncluded(modelVersion, lognormalMassAction)) {
    normp=W_UNIFORM();
    Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
    expectedNumInfections=Infectionrate;
        //Bad (duplicated) code
  } else if( isOptionIncluded(modelVersion, lognormalMassActionPlusPreImm)) {
    normp=W_UNIFORM();
    Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
    survivalOfInoculum=(1.0+pow((human.getCumulativeEIRa()/Xstar_p), gamma_p));
    survivalOfInoculum=Simm+(1.0-Simm)/survivalOfInoculum;
    survivalOfInoculum=survivalOfInoculum*(Sinf+(1-Sinf)/(1+effectiveEIR/Estar));
    expectedNumInfections=survivalOfInoculum*Infectionrate;
  } else {
    //The default model is that in Smith et al, AJTMH 2006 75 Suppl 2
    survivalOfInoculum=(1.0+pow((human.getCumulativeEIRa()/Xstar_p), gamma_p));
    survivalOfInoculum=Simm+(1.0-Simm)/survivalOfInoculum;
    survivalOfInoculum=survivalOfInoculum*(Sinf+(1-Sinf)/(1+effectiveEIR/Estar));
    expectedNumInfections=survivalOfInoculum*effectiveEIR*interval*baseAvailToMos;
  }
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if ( isOptionIncluded(vaccineType, preerythrocytic_reduces_h)) {
    expectedNumInfections=expectedNumInfections*(1-efficacy);
  }
  //Update pre-erythrocytic immunity
  if (isOptionIncluded(modelVersion,transHet) || isOptionIncluded(modelVersion,comorbTransHet) || isOptionIncluded(modelVersion,transTreatHet) || isOptionIncluded(modelVersion, tripleHet)) {
    human.updateCumulativeEIRa((double)(interval*effectiveEIR*baseAvailToMos));
  } else {
    human.updateCumulativeEIRa((double)interval*effectiveEIR);
  }
  return expectedNumInfections;
}

double TransmissionModel::getExpectedNumberOfInfections(double efficacy, double baseAvailToMos, double cumulativeEIRa,  double age_adj_EIR) {
  double normp;
  double survivalOfInoculum;
  double ExpectedInfectionRate;
  double Infectionrate;
  //The age-adjusted EIR, possibly adjusted for bed nets.
  double expectedNumInfections;
  
  ExpectedInfectionRate=age_adj_EIR*baseAvailToMos*susceptibility*interval;
  if ( isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
    Infectionrate=(W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
    expectedNumInfections=Infectionrate;
  }
  else if( isOptionIncluded(modelVersion, lognormalMassAction)) {
    normp=W_UNIFORM();
        Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
        expectedNumInfections=Infectionrate;
        //Bad (duplicated) code
  }
  else if( isOptionIncluded(modelVersion, lognormalMassActionPlusPreImm)) {
    normp=W_UNIFORM();
    Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
    survivalOfInoculum=(1.0+pow((cumulativeEIRa/Xstar_p), gamma_p));
    survivalOfInoculum=Simm+(1.0-Simm)/survivalOfInoculum;
    survivalOfInoculum=survivalOfInoculum*(Sinf+(1-Sinf)/(1+age_adj_EIR/Estar));
    expectedNumInfections=survivalOfInoculum*Infectionrate;
  }
  else {
    //The default model is that in Smith et al, AJTMH 2006 75 Suppl 2
    survivalOfInoculum=(1.0+pow((cumulativeEIRa/Xstar_p), gamma_p));
    survivalOfInoculum=Simm+(1.0-Simm)/survivalOfInoculum;
        survivalOfInoculum=survivalOfInoculum*(Sinf+(1-Sinf)/(1+age_adj_EIR/Estar));
        expectedNumInfections=survivalOfInoculum*age_adj_EIR*interval*baseAvailToMos;
  }
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if ( isOptionIncluded(vaccineType, preerythrocytic_reduces_h)) {
    expectedNumInfections=expectedNumInfections*(1-efficacy);
  }
  return expectedNumInfections;
}

double TransmissionModel::calculateEIR(int simulationTime){
  // Calculates EIR (in adults), based on vectorial capacity or looks up EIR in
  // the input data.
  // time: Time since start of simulation .

  // where the full model, with estimates of human mosquito transmission is in use, use this:
  switch (simulationMode) {
    case equilibriumMode:
      return EIR[modIntervalsPerYear(simulationTime) - 1];
      break;
    case transientEIRknown:
      // where the EIR for the intervention phase is known, obtain this from
      // the intEIR array
      return intEIR[Simulation::timeStep - 1];
      break;
    case dynamicEIR:
      if ( Simulation::timeStep ==  1) {
        return EIR[modIntervalsPerYear(simulationTime) - 1];
      } else if ( ITN) {
        double Pc0=pow(Pu0, c);
        double Puz=Pu0-z*(Pu0-Pu1);
        double Pcz=pow(Puz, c);
        double dz=(1.0-Pu0)/(1.0-Puz);
        double s0t=initialKappa[modIntervalsPerYear(simulationTime-nearbyint(nspore)) - 1] * Pc0/(1-Pu0);
        double szt=kappa[modIntervalsPerYear(simulationTime-nearbyint(nspore)) - 1] * Pcz/(1-Puz);
        return EIR[modIntervalsPerYear(simulationTime) - 1] *
            dz * szt / s0t;
      } else {
        return EIR[modIntervalsPerYear(simulationTime) - 1] *
            kappa[modIntervalsPerYear(simulationTime-nearbyint(nspore)) - 1] /
            initialKappa[modIntervalsPerYear(simulationTime-nearbyint(nspore)) - 1];
      }
      break;
    default:
      cerr << "Invalid simulation mode" << endl;
      throw 0;	// Anything else.. don't continue silently
  }
}


void TransmissionModel::logDFTThreeModeSmooth (double* smoothArray,
                                               double* originalArray,
                                               int SALength,
                                               int OALength) {
    /*
    ***************************************************************
    **************** logDFTThreeModeSmoothExpand ******************
    ***************************************************************
    *  Given a positive array, originalArray, of length OALength,
    *  this routine exponentiates the inverse discrete Fourier 
    * tranform of the first three modes of the natural logarithm of 
    * the array to smooth out the array to produce smoothArray of 
    * length SALength.
    *
    * All elements of originalArray are assumed to be strictly
    * positive.
    *
    * smoothArray is an OUT parameter.
    * originalArray, SALength and OALength are IN parameters.
    */

    // Frequency
    double foa = 1.0/OALength;
    double woa = 2*M_PI * foa;
    double wsa = 2*M_PI/SALength;
    
    double tempsuma0 = 0.0;
    double tempsuma1 = 0.0;
    double tempsumb1 = 0.0;
    double tempsuma2 = 0.0;
    double tempsumb2 = 0.0;

    // Calculate first three Fourier modes
    for (int t=0; t < OALength; t++) {
        double yt = log(originalArray[t]);
        double woa_t = woa*t;
        tempsuma0 = tempsuma0+yt;
        tempsuma1 = tempsuma1 + (yt*cos(woa_t));
        tempsumb1 = tempsumb1 + (yt*sin(woa_t));
        tempsuma2 = tempsuma2 + (yt*cos(2*woa_t));
        tempsumb2 = tempsumb2 + (yt*sin(2*woa_t));       
    }
    
    // Fourier Coefficients
    double a0 = (  foa)*tempsuma0;
    double a1 = (2*foa)*tempsuma1;
    double b1 = (2*foa)*tempsumb1;
    double a2 = (2*foa)*tempsuma2;
    double b2 = (2*foa)*tempsumb2;
    
    // Calculate inverse discrete Fourier transform
    for (int t=0; t < SALength; t++){
        double wsa_t = wsa*(t+1);
        smoothArray[t] = 
            exp(a0 + a1*cos(wsa_t) + b1*sin(wsa_t) + 
                a2*cos(2*wsa_t) + b2*sin(2*wsa_t));
    }

#   ifdef TransmissionModel_PrintSmoothArray
        PrintArray(fnametestentopar, "SmoothArray", smoothArray, SALength);
#   endif
}



void TransmissionModel::calcInverseDFTExp(double* tArray, int aL, double* FC, int FCL) {
    /*
    ***************************************************************
    ********************** calcInverseDFTExp **********************
    ***************************************************************
    *  Given a sequence of Fourier coefficients, FC, of length, FCL,
	*  this routine calculates the exponent of the inverse discrete
	*  Fourier transform into an array, Tarray, of length, aL.
    *
	*  Note that FCL is assumed to be an odd number.
	*  
	* tArray is an OUT parameter.
	* aL, FC, and FCL are IN parameters.
    */

    // Period
    double P = (double)aL;
    // Frequency
    double w = 2*M_PI/P;
    
    if((FCL%2)==0){
        printf("The number of Fourier coefficents should be odd.\n");
        getchar();
    }

    // Number of Fourier Modes.
    int Fn = (FCL-1)/2;

    // Calculate inverse discrete Fourier transform
    for (int t=1; t<=aL; t++){
       double temp = FC[0];
        if(Fn>0){
            for(int n=1;n<=Fn;n++){
                temp = temp + FC[2*n-1]*cos(n*w*t) + FC[2*n]*sin(n*w*t);
            }
        }
        tArray[t-1] = exp(temp);
    }
}

void TransmissionModel::rotateArray(double* rArray, int aLength, double rAngle) {
    /*
    ***************************************************************
    ************************ rotateArray **************************
    ***************************************************************
    *  Given an array, rArray, of length aLength, the routine rotates
    *  the array clockwise by rAngle.
    *
    * rArray is an IN/OUT parameter.
    * aLength and rAngle are IN parameters.
    */

    double* tempArray = (double *)malloc((aLength)*sizeof(double));

    int rotIndex = (int)((rAngle*aLength)/(2*M_PI));

#   ifdef TransmissionModel_PrintRotateArray
        PrintArray(fnametestentopar, "PrerotationArray", rArray, aLength);
#   endif


    for (int i=0;i<aLength; i++) {
        int tIndex = (i+rotIndex) % aLength;
        tempArray[tIndex] = rArray[i];
    }
    
    for (int i=0; i<aLength; i++){
        rArray[i] = tempArray[i];
    }

#   ifdef TransmissionModel_PrintRotateArray
        PrintArray(fnametestentopar, "PostrotationArray", rArray, aLength);
#   endif
    free(tempArray);
}


void TransmissionModel::PrintArray(char fntestentopar[], char vectorname[], double* v, int n){
/********************************************************************/
/* PrintArray() prints the given (C) array to the given file.
 * 
 * The array, v, of doubles is assumed to be of length n.
 * All parameters are IN parameters.
 */

  FILE* fpp = fopen(fntestentopar, "a");

  for (int i=0; i < n; i++){
    fprintf(fpp, "%s(%d) = %f; \n", vectorname, i+1, v[i]);
  }
  fclose(fpp);
}

