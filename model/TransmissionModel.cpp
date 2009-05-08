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
//This is needed to get symbols like M_PI with MSVC:
#define _USE_MATH_DEFINES

#include "TransmissionModel.h"
#include <math.h> 
#include "global.h" 
#include <gsl/gsl_vector.h> 
#include "intervention.h" 
#include <iostream>
#include "inputData.h"
#include "TransmissionModel/NonVector.h"
#include "TransmissionModel/Vector.h"
#include "TransmissionModel/PerHost.h"
#include "simulation.h"
#include "summary.h"

//static (class) variables
const double TransmissionModel::agemin[nwtgrps] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60 };
const double TransmissionModel::agemax[nwtgrps] = { 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 19.99, 24.99, 29.99, 39.99, 49.99, 59.99, 60.99 };
const double TransmissionModel::bsa_prop[nwtgrps] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

TransmissionModel* TransmissionModel::createTransmissionModel () {
  // EntoData contains either a list of at least one anopheles or a list of at
  // least one EIRDaily.
  if (getEntoData().getAnopheles().size())
    return new VectorTransmission();
  else
    return new NonVectorTransmission();
}

TransmissionModel::TransmissionModel(){
  strcpy(fnametestentopar,"output_ento_para.txt");
  
  EIRRotateAngle = M_PI_2;
  EIPDuration = get_EipDuration();
  PerHostTransmission::BaselineAvailabilityShapeParam=getParameter(Params::BASELINE_AVAILABILITY_SHAPE);	// also set in Human
  
  kappa = (double *) malloc(((Global::intervalsPerYear))*sizeof(double));
  initialKappa = (double *) malloc(((Global::intervalsPerYear))*sizeof(double));
  EIR = new double[Global::intervalsPerYear];
  origEIR = new double[Global::intervalsPerYear];
  
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

  initAgeExposureConversion();
}

TransmissionModel::~TransmissionModel () {
  free(kappa);
  free(initialKappa);
  delete [] EIR;
  delete[] origEIR;
  free(FCEIR);
}

void TransmissionModel::updateKappa (double sumWeight, double sumWt_kappa) {
  int tmod = (Simulation::simulationTime - 1) % Global::intervalsPerYear;
  //Prevent NaNs
  if (sumWeight == 0.0) {
    kappa[tmod] = 0.0;
    cerr << "sW.eq.0" << endl;
  }
  else {
    kappa[tmod] = sumWt_kappa / sumWeight;
  }
  
  //Calculate time-weighted average of kappa
  if (tmod == 0) {
    _sumAnnualKappa = 0.0;
  }
  _sumAnnualKappa += kappa[tmod] * Global::interval * EIR[tmod];
  if (tmod + 1 == Global::intervalsPerYear) {
    if (annualEIR == 0) {
      _annualAverageKappa=0;
      cerr << "aE.eq.0" << endl;
    }
    else {
      _annualAverageKappa = _sumAnnualKappa / annualEIR;
    }
  }
}

void TransmissionModel::summarize (Summary& summary) {
  summary.setNumTransmittingHosts(kappa[(Simulation::simulationTime-1) % Global::intervalsPerYear]);
  summary.setAnnualAverageKappa(_annualAverageKappa);
}

double TransmissionModel::getRelativeAvailability(double ageyrs) {
  return ageSpecificRelativeAvailability[getAgeGroup(ageyrs)];
}

void TransmissionModel::initAgeExposureConversion(){
  for (size_t i=0; i<nages; i++) {
    ageSpecificRelativeAvailability[i] = bsa_prop[i] / (1-bsa_prop[i]);
  }

  /* Not used.
  // The following code calculates average bites for children <6 years old
  // relative to adults. This is for analysis of Saradidi (Beier et al) data.

  // avbites_6 is the exposure of children < 6 years old relative to adults
  double avbites_6=0.0;
  
  for (int i = 0; agemin[i] <  6; ++i) {
    if (agemax[i] <= 0.5)	// NOTE: never true − is it obsolete?
      continue;
    
    //agemin0 is the lower bound of the youngest age group used in this calculation
    double agemin0=agemin[i];
    if (agemin0 < 0.5)
      agemin0=0.5;
    
    avbites_6 += ageSpecificRelativeAvailability[i] * (agemax[i]-agemin0)/5.5;
  }
  biteratio_6 = avbites_6 / (1.0-avbites_6);
  */
}

void TransmissionModel::logDFTThreeModeSmooth (double* smoothArray,
                                               double* originalArray,
                                               int SALength,
                                               int OALength)
{
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
    // Period
  double P = (double)aL;
    // Frequency
  double w = 2*M_PI/P;
    
  if((FCL%2)==0){
    //NOTE: should throw xml_scenario_error if/when FCEIR is moved to the XML
    throw logic_error("The number of Fourier coefficents should be odd.");
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
  double* tempArray = (double *)malloc((aLength)*sizeof(double));

  int rotIndex = (int)((rAngle*aLength)/(2*M_PI));

#   ifdef TransmissionModel_PrintRotateArray
  PrintArray(fnametestentopar, (char*)"PrerotationArray", rArray, aLength);
#   endif


  for (int i=0;i<aLength; i++) {
    int tIndex = (i+rotIndex) % aLength;
    tempArray[tIndex] = rArray[i];
  }
    
  for (int i=0; i<aLength; i++){
    rArray[i] = tempArray[i];
  }

#   ifdef TransmissionModel_PrintRotateArray
  PrintArray(fnametestentopar, (char*)"PostrotationArray", rArray, aLength);
#   endif
  free(tempArray);
}


void TransmissionModel::PrintArray(char fntestentopar[], char vectorname[], double* v, int n){
  FILE* fpp = fopen(fntestentopar, "a");

  for (int i=0; i < n; i++){
    fprintf(fpp, "%s(%d) = %f; \n", vectorname, i+1, v[i]);
  }
  fclose(fpp);
}

size_t TransmissionModel::getAgeGroup (double age) {
  for (size_t i = 0; i < nages; ++i) {
    if (agemax[i] > age)
      return i;
  }
  return nages-1;	// final category
}


// -----  checkpointing  -----

void TransmissionModel::write(ostream& out) const {
  out << annualEIR << endl;
  for (int i = 0; i < Global::intervalsPerYear; ++i)
    out << kappa[i] << endl;
  out << _annualAverageKappa << endl;
  out << _sumAnnualKappa << endl;
}
void TransmissionModel::read(istream& in) {
  in >> annualEIR;
  for (int i = 0; i < Global::intervalsPerYear; ++i)
    in >> kappa[i];
  in >> _annualAverageKappa;
  in >> _sumAnnualKappa;
}
