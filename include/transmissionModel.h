#ifndef Hmod_TransmissionModel 
#define Hmod_TransmissionModel 
#define _USE_MATH_DEFINES 

#include "human.h"
#include "string.h"

// Define these to print out various arrays:
//#define TransmissionModel_PrintOrigEIR
//#define TransmissionModel_PrintEIRaIDFT
//#define TransmissionModel_PrintSmoothArray
#define TransmissionModel_PrintRotateArray

 
//! Abstract base class, defines behaviour of transmission models
class TransmissionModel{ 

public: 

 TransmissionModel();
 virtual ~TransmissionModel() {}; 
 
  //! Reads all entomological parameters from the input datafile. 
 void initEntoParameters();

 //! get the number of infections for a specific human at a given time step 
 /*! 
 1. Calculates h from the EIR measured on adults where 
 h is the expected number of epidemiological inoculations 
 2 Calculates the updated values of the pre-erythrocytic exposure and 
 passes this back to the calling routine 
 Requires the five-day EIR, adjusted for age as input. 
 cumEIR: is the pre-erythrocytic exposure; 
 efficacy: Efficacy of a pre-erythrocytic vaccine 
 \param *cumEIR cumulative EIR (measures pre-erthrocytic immunity) 
	\param efficacy efficacy of a pre-erythrocytic vaccine 
	\param age_adj_EIR Expected number of inoculations adjusted for age of the host 
	\param baseAvailToMos Host-specific availability 
 */ 
 virtual double getExpectedNumberOfInfections (double efficacy, double baseAvailToMos, double cumulativeEIRa,  double age_adj_EIR);
 virtual double getExpectedNumberOfInfections (Human& human, double age_adj_EIR);

//! Calculates the adjustment for body size in exposure to mosquitoes 
/*! 
The bites are assumed proportional to average surface area for hosts of the given age. 
Linear interpolation is used to calculate this from the input array of surface areas. 
\param ageyrs age in years 
\return the ratio of bites received by the host to the average for an adult 
	 */ 
double getRelativeAvailability (double ageyrs); 


double getWeight(double age);

//! initialise the main simulation 
virtual void initMainSimulation (int populationSize)=0; 
 
//!Deallocate memory for TransmissionModel parameters and clean up (TODO: change to destructor) 
void clearTransmissionModelParameters (); 

//kappa[] is the probability of infection of a mosquito at each bite 
    double *kappa; 

 //initialKappa[] is the value of kappa during the pre-intervention phase  
    double *initialKappa; 

//TODO: at least the following should be made private
//protected: 
//functions used by the constructor 

//! read in the EIR from the .xml 
/*! should have different versions for the VectorControl and NoVectorControl classes 
*/
void inputEIR (); 
 
//! EIR per time step during the pre-intervention phase 
double *EIR; 
 
//FUNCTIONS THAT SHOULD BE USED BY getExpectedNumberOfInfections (renamed as getNinfections ???) 

//! get the EIR at a given time step (getPopulationInoculationRate ???) 
double calculateEIR(int simulationTime); 

// There should be a function declared here that does the core gets

//Variables used by calculateEIR 
//! total annual EIR 
double annualEIR; 
 
protected:
//! The duration of sporogony in time steps 
double nspore; 

private:
/* 
    nwtgrps is the number of age groups for which expected weights are read in 
    for use in the age adjustment of the EIR 
*/
static const int nwtgrps= 27; 

//! average number of bites for the age as a proportion of the maximum 
double ageSpecificRelativeAvailability[nwtgrps]; 
//! Relative weights by age group
 /* 
Relative weights, based on data in InputTables\wt_bites.csv 
The data are for Kilombero, Tanzania, taken from the Keiser et al (diploma thesis). The original source 
was anthropometric studies by Inez Azevedo 
Reads in weights by age group. The weights are expressed as proportions of 0.5*those in 
the reference age group 
 */ 
static const double wtprop[nwtgrps]; 

//! Cutpoints of the age categories (minima) used for storing relative weights? surface areas?
static const double agemin[nwtgrps]; 

//! Cutpoints of the age categories (maxima) used for storing relative weights? surface areas?
static const double agemax[nwtgrps];  

//! Proportionate body surface area
 /* 
The body surface area is expressed as proportions of 0.5*those in 
the reference age group.In some models we have used calculations of weight and in others surface area, based on 
Mosteller RD: Simplified Calculation of Body Surface Area. N Engl J Med 1987 Oct 22;317(17):1098 (letter) 
These values are retained here should they be required for future comparisons 
 */ 
static const double bsa_prop[nwtgrps]; 

// ratio of the number of bites received relative to the number received at age 6 years 
double biteratio_6; 

protected:
static const int ifrotateEIR = 0;	// TODO: Move to XML.
							// Flag to rotate EIR by a given number of days
							// to check the calculations for Kappa.
static const int ifUseFC = 0;		// TODO: Move to XML.
							// Flag to use Fourier coefficients to create EIR (instead of time series data).
							// Right now we do not link this to FTSmoothEIR - but these definitions should
							// be linked.
   double* origEIR;	// Original EIR - if we smooth the EIR using the Fourier transform.
   double* FCEIR;	// Fourier coefficients for the EIR.
   int	FCEIRX;		// Number of Fourier coefficients used to calculate EIR.
   double EIRRotateAngle;	// Angle to rotate EIR: Should be between 0 and 2Pi.
   char fnametestentopar[30];
   static const int FTSmoothEIR = 0;	// TODO: Move to XML: 1 to smooth EIR using an approximate DFT
						//                    0 to do nothing.

void logDFTThreeModeSmooth (double* smoothArray, double* originalArray, int SALength, int OALength); 
void calcInverseDFTExp(double* tArray, int aL, double* FC, int FCL); 
void rotateArray(double* rArray, int aLength, double rAngle); 
void PrintArray(char fntestentopar[], char vectorname[], double* v, int n);


private: 
//! multiplier used to calculate a positive EIR value where the measured value is zero
/* 
    0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when it is too low 
    to be measureable. The value of 0.01 was old pv(30) Now a constant. 
    0.01 was old pv(30) Now a constant. 
*/ 
    static const double min_EIR_mult; 

//! Number of days contributing to each EIR estimate for pre-intervention 
    int *no; 

//VARIABLES INCLUDED IN CORE GETs of number of infections 
//! The average proportion of bites from sporozoite positive mosquitoes resulting in infection. 
 /*! 
 This is computed as 0.19 (the value S from a neg bin mass action model fitted 
	to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
 0.5-6 year old child (as per Saradidi) to adult) 
 \sa getExpectedNumberOfInfections() 
*/ 
static const double susceptibility; 

//!Steepness of relationship between success of inoculation and Xp in Phase A model 
 /*! 
	\sa getExpectedNumberOfInfections(),Sinf,Simm,Xstar_p,Estar 
*/ 
double gamma_p; 
 
//!Lower limit of success probability of inoculations at high exposure in Phase A model 
 /*! 
	\sa getExpectedNumberOfInfections(),gamma_p,Simm,Xstar_p,Estar 
*/ 
 double Sinf; 
 
//!Lower limit of success probability of inoculations in immune individuals in Phase A model 
 /*! 
	\sa getExpectedNumberOfInfections(),gamma_p,Sinf,Xstar_p,Estar 
*/ 
 double Simm; 
 
//!Critical value of cumulative number of entomologic inoculations in Phase A model 
 /*! 
	\sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Estar 
*/ 
 double Xstar_p; 
 
//!Critical value of EIR in Phase A pre-erythrocytic model 
 /*! 
	\sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Xstar_p 
*/ 
 double Estar; 
 

//! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr. 
double InfectionrateShapeParam; 

//! Variance of Infection Rate according to fielddata 
static const double totalInfectionrateVariance; 

//! Shape parameter of the distribution of the baseline availability
/* 
Same BaselineAvailabilityShapeParam as in human.f, but we have to redefine it here, since mod_human 
uses MOD_entomology 
*/ 
 double BaselineAvailabilityShapeParam; 

//functions used by initMainSimulation
 
//! Number of age groups for which the surface area calculations apply 
static const int nages= 22;


/*
We add entomological model parameters here. 
MosqEmergeRate is the emergence rate of mosquitoes over one period. 
This is described in more detail below. 
*/ 
double mosqEmergeRate[daysInYear]; 

 
// we expect to put the following into NoVectorControl class 
void updateEIR (int day, double EIRdaily); 
double averageEIR (); 

//TODO: the entire code for specifying availability should be part of the human
//! initialisation of the vector of expected surface area as a function of age
void initAgeExposureConversion(); 

//! The maximum number of daily EIR values specified for intervention phase 
 /*! 
	The EIR is input on a daily basis even when the interval length is >1 day 
	Average EIR per interval is computed in initEntoParameters 
	\sa initEntoParameters 
 */ 
static const int maxDurIntPhaseEIR= 1500; 
 
//! The maximum number of intervals in the intervention phase 
int maxIntervals; 

//ino() Number of days contributing to each EIR estimate (post intervention) 
int *ino; 

int inoX; 

 //intEIR() EIR per time interval during the intervention period 
double *intEIR; 

int intEIRX; 

}; 
 
#endif 
