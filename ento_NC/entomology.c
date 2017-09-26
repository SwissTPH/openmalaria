#include "entomology.h"
#include "bridge.h"
/*

 This file is part of mcdnsa.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute
 
 mcdnsa is free software; you can redistribute it and/or modify
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
//maxDurIntPhaseEIR is the maximum number of days for which an intervention phase EIR is specified	
	int FTSmoothEIR = 0;	// TODO: Move to XML: 1 to smooth EIR using an approximate DFT
						//                    0 to do nothing.
   	int ifrotateEIR = 0;	// TODO: Move to XML.
							// Flag to rotate EIR by a given number of days
							// to check the calculations for Kappa.
	int ifUseFC = 0;		// TODO: Move to XML.
							// Flag to use Fourier coefficients to create EIR (instead of time series data).
							// Right now we do not link this to FTSmoothEIR - but these definitions should
							// be linked.
	double EIRRotateAngle = M_PI/2;	// Angle to rotate EIR: Should be between 0 and 2Pi.

	// File name where we send output of entomological model.
	char fnametestentopar[30] = "output_ento_para.txt";	
   double nspore;
   int maxIntervals;
   double gamma_p;
   double Sinf;
   double Simm;
   double Xstar_p;
   double Estar;
   double *kappa;
     int kappaX;
   double *initialKappa;
     int initialKappaX;
   double *EIR;
    int EIRX;
   double* origEIR;	// Original EIR - if we smooth the EIR using the Fourier transform.
   double* FCEIR;	// Fourier coefficients for the EIR.
   int	FCEIRX;		// Number of Fourier coefficients used to calculate EIR.
   double annualEIR;
   int *no;
     int noX;
   int *ino;
     int inoX;
   double *intEIR;
     int intEIRX;
   double btprop[nwtgrps];
   double biteratio_6;
   double InfectionrateShapeParam;
   double BaselineAvailabilityShapeParam;
    double agemin[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60 };
    double agemax[] = { 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 19.99, 24.99, 29.99, 39.99, 49.99, 59.99, 60.99 };
    double wtprop[] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
    double bsa_prop[] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
   double mosqEmergeRate[daysInYear];
void initEntoParameters () {
/*

Read all entomological parameters from the input datafile. 

*/

//Product of measured EIR, susceptibility and Time inerval
    double gsi;
//constant defining the constraint for the Gamma shape parameters
    double r_square_Gamma;
//constant defining the constraint for the LogNormal parameters
    double r_square_LogNormal;
    nspore=get_nspore();
    gamma_p=get_parameter(5);
    Sinf=1-exp(-get_parameter(1));
    Simm=get_parameter(3);
    Estar=get_parameter(2);
    Xstar_p=get_parameter(4);
    BaselineAvailabilityShapeParam=get_parameter(16);
    maxIntervals=maxDurIntPhaseEIR/interval;
    no = malloc(((intervalsPerYear))*sizeof(int));
    noX = intervalsPerYear;
    kappa = malloc(((intervalsPerYear))*sizeof(double));
    kappaX = intervalsPerYear;
    initialKappa = malloc(((intervalsPerYear))*sizeof(double));
    initialKappaX = intervalsPerYear;
    EIR = malloc(((intervalsPerYear))*sizeof(double));
    EIRX = intervalsPerYear;
	origEIR = malloc(((intervalsPerYear))*sizeof(double));
    ino = malloc(((maxIntervals))*sizeof(int));
    inoX = maxIntervals;
    intEIR = malloc(((maxIntervals))*sizeof(double));
    intEIRX = maxIntervals;
    gsi=1.0;
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

    /*
    r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
    r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
    */
    r_square_Gamma=0.649;
    //such that r_square_LogNormal =0.5
    r_square_LogNormal=log(1.0+r_square_Gamma);
    //TODO: Sanity check for sqrt and division by zero
    if ( isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
        InfectionrateShapeParam=(BaselineAvailabilityShapeParam+1.0)/(r_square_Gamma*BaselineAvailabilityShapeParam-1.0);
        InfectionrateShapeParam=max(InfectionrateShapeParam, 0.0);
    }
    else if( isOptionIncluded(modelVersion, lognormalMassAction) ||  isOptionIncluded(modelVersion, lognormalMassActionPlusPreImm)) {
        InfectionrateShapeParam=sqrt(r_square_LogNormal-1.86*pow(BaselineAvailabilityShapeParam, 2));
        InfectionrateShapeParam=max(InfectionrateShapeParam, 0.0);
    }
    inputEIR();
    surfaceAreaAgeConversion();
}

void calMosqEmergeRate () {
    /*
    ***************************************************************!
    ********************** calMosqEmergeRate ********************!!
    
     Add subroutine here that calls C function that initializes 
     the system. We can call this subroutine calMosqEmergeRate
     This routine passes the basic entomological parameters (that 
     we assume will already have been read, the EIR, and the human 
     infectivity to mosquitoes (all for one type of host) to a C 
     function that then calculates the mosquito emergence rate 
     (over a one year period). This rate will then be used in the 
     main simulation.
    
     For now, we just call testFortranCCommunication to see how
     well we can pass arrays and matrices between the two.
    */

	/*
     Number of types of hosts. 
     Dimensionless.
     $n$ in model. Scalar.
     Is equal to 1 in initalization.
    */
    int nHostTypesInit = 1;

    /*
     Number of type of malaria-susceptible hosts. 
     Dimensionless.
     $m$ in model. Scalar.
     Is equal to 1 in initalization.
    */
    int nMalHostTypesInit = 1;
//;
    /*
     Number of hosts of each type.
     Units: Animals. 
     $N_i$ in model. Matrix of size $n \times \theta_p$.
     We assume that the size of the one group in initialization is
     fixed over the cycle.
     Mathematically, we require this parameter to be a positive
     real number, so although this will typically be a natural 
     number, it is not restricted to being one.
    */
    double popSizeInit;
    /*
     Infectivity of hosts to mosquitoes.
     Dimensionless.
     $K_vi$ in model. Matrix of size $n \times \theta_p$.
     In initialization, there is only one time of host.
     This is taken directly from initialKappa(i) from
     entomology.f
     Also needs to be converted to real*8 
     (and to the correct length).
    */
	double humanInfectivityInit[daysInYear];
    /*
     Duration of the extrinsic incubation period.
     Units: Days.
     $\theta_s$ in model. Scalar.
     We calculate this from nspore and interval.
    */
    int EIPDuration;
    /*
    ***************************************************************!
    ***************************************************************!
    ** We enter parameters here that will later be moved to .xml **!
    ***************************************************************!
     Note that these parameters are only for one group of humans.
     We need to get parameters for multiple groups later.
     Some parameters - that should be vectors of length daysInYear
     we will enter as scalars and assume they are fixed over the 
     year. Although, in the theoretical model, we assume that they
     may vary periodically over one year, we code them as though 
     they are constant. We can change that later if they need to 
     vary.
     We assume that we will read them from the xml file as doubles
     so that they will be passed into Fortran as real*8. We 
     therefore define them here as real*8. The results, then, when
     used in the main simulation can be converted to reals.
     When they are defined in the xml file, we can rearrange 
     parameters here so that one group is taken directly from the
     xml file and the other group from other Fortran routines.
     And we still initialize the parameters here, but do not
     define them.
     Duration of resting period for mosquito.
     This will be defined in the xml file as a real.
     We will define it as a integer - and convert after multiplying
     by the length of the interval. The unit is days.
     $\tau$ in model. Scalar.
    */
    int mosqRestDuration;
    /*
     Availability rate of hosts to mosquitoes.
     Units: 1/(Animals * Time)
     $\alpha_i$ in model. Matrix of size $n \times \theta_p$.
     In initialization, there is only one time of host.
     We also assume that this does not change over the cycle.
     It appears that this is not defined in Fortran. We need
     to define this quantity ourselves.
     This adds a further degree of freedom - the value should be 
     chosesn carefully as there will be an inverse relationhip
     between $\sum_{i=1}^n \alpha_i N_i$ and $N_{v0}$.
    */
    double hostAvailabilityRateInit;
    /*
     Mosquito death rate while host-seeking. The unit is 1/Days.
     $\mu_{vA}$ in model. Vector of length $\theta_p$.
    */
    double mosqSeekingDeathRate;
    /*
     Duration of host-seeking per day. This is the maximum fraction 
     of a day that a mosquito would spend seeking if it is
     unsuccessful in finding a host. Dimensionless
     This parameter is defined in the model as a
    */
    double mosqSeekingDuration;
    /*
     Probability of a mosquito biting a host given that it has 
     encountered the host. 
     Dimensionless.
     $P_{B_i}$ in model. Matrix of size $n \times \theta_p$.
     For now, we assume this does not change over the cycle.
    */
    double mosqProbBiting;
    /*
     Probability of a mosquito finding a resting side given that it
     has bitten a host.
     Dimensionless.
     $P_{B_i}$ in model. Matrix of size $n \times \theta_p$.
     For now, we assume this does not change over the cycle.
    */
    double mosqProbFindRestSite;
    /*
     Probability of a mosquito surving the resting period given that
     it has found a resting site.
     Dimensionless.
     $P_{D_i}$ in model. Matrix of size $n \times \theta_p$.
     For now, we assume this does not change over the cycle.
    */
    double mosqProbResting;
    /*
     Probability of a mosquito ovipositing and returning to
     host-seeking given that it has survived the resting period.
     Dimensionless.
     $P_{E_i}$ in model. Matrix of size $n \times \theta_p$.
     For now, we assume this does not change over the cycle.
    */
    double mosqProbOvipositing;
    /*
    ***************************************************************! 
    ***************************************************************!
    *********** Other main input and output parameters ************!
    ***************************************************************!
     The mosquito emergence rate. 
     Units: Mosquitoes/Time
     $N_{v0}$ in model. Vector of length $\theta_p$.
     This subroutine calls a function in C to calculate this rate.
     It is defined globally as mosqEmergeRate.
     The initial estimate of the mosquito emergence rate. 
     Units: Mosquitoes/Time
     Not explicitly included in model. Vector of length $\theta_p$.
     This initial estimate is used by a root finding algorithm
     to calculate the mosquito emergence rate.
    */
    double mosqEmergeRateInitEstimate[daysInYear];
    /*
     The entomological inoculation rate.
     Units: 1/Time
     $\Xi_i$ in model. Matrix of size $n \times \theta_p$.
     We use this to calculate the mosquito emergence rate.
     During the initialization, it a vector with the length of the 
     annual period.
    */
	double EIRInit[daysInYear];

	int ifUseNv0Guess = 0;	// Use a predefined array for the initial
							// mosquito emergence rate.
							// Perhaps calculated in a different
							// iteration.

    /*
    ***************************************************************!
    ***************************************************************!
    *************** Other parameters and variables ****************!
    ***************************************************************!
    */
    int i;
    int printvectors;
	int din = daysInYear;
    int hti = nHostTypesInit;
    int mhti = nMalHostTypesInit;

	int smoothfullyearEIR = 1;	// If smoothing the EIR, smooth it over the
								// entire year (and not in increments of 
								// intervalLength.

	int ifprintKvi = 0;	// Print the initial human infectivity to mosquitoes.
	char shortKviname[15] = "ShortKvi";
	char longKviname[15] = "LongKvi";

	int ifprintEIR = 0;	// Print the initial EIR.
	char origEIRname[15] = "OrigEIR";
	char shortEIRname[15] = "ShortEIR";
	char longEIRname[15] = "LongEIR";
    /*
    ***************************************************************!
    ***************************************************************!
    *************** Begin code of subroutine here *****************!
    ***************************************************************!
     Be careful between real and real*8.
     Set values of parameters that we define here
     (that will later be moved to the xml file).
    */
    mosqRestDuration=3;			// TODO: Move to XML
    mosqSeekingDeathRate=1.6;	// TODO: Move to XML
    mosqSeekingDuration=0.33;	// TODO: Move to XML
    mosqProbBiting=0.95;		// TODO: Move to XML
    mosqProbFindRestSite=0.95;	// TODO: Move to XML
    mosqProbResting=0.94;		// TODO: Move to XML
    mosqProbOvipositing=0.93;	// TODO: Move to XML
    // We set the host availability relative to the population size.
    hostAvailabilityRateInit=7.0/(npeople);	// Move absolute availability to XML
    // Set scalars currently read from xml/other Fortran routines.
    popSizeInit=(npeople);
    EIPDuration=floor(nspore*interval);

	


    /*
     Now we have to deal with 
     - humanInfectivityInit - we take from initialKappa
     - EIRinit - we take from EIR.
    ************ Make sure that EIR is the correct one. ***********!
     I'm not sure how we should check this - but we should look into
     it more carefully at some point in time. For now, we can 
     continue with this.
     We first create arrays of length intervalsperYear for all three.
     We then convert them to length daysInYear.
     Save the human infectivity to mosquitoes from simulation of one
     lifetime to initialKappa. We will then use this to create 
     humanInfectivityInit (of size daysInYear).
     Do not mess with initialKappa - this may be used elsewhere.
    */
    for ( i=1;i<=intervalsPerYear; i++) {
        initialKappa[i - 1]=kappa[i - 1];
    }
    /*
     We need to convert arrays of length intervalsPerYear to length
     daysInYear.
    */

	// We need to decide how we deal with the EIR - if we smooth it out
	// over the entire year, or leave it constant over the interval length.
	// I think it would be better ot smooth it over the full year. 
	// It may not be fully accurate - but we are in any case going to lose
	// some accuracy over the difference between the time step of the human
	// simulation model and the mosquito transmission model.
	// Note that we smooth over the entire year, we are slightly shifting
	// the EIR a little bit to the right.
	// We have a lot of if statements and flags here. We need to clean this 
	// up eventually when the flags are moved to XML.
	if(smoothfullyearEIR==1){
		if(FTSmoothEIR==1){
			printf("Smoothing and expanding EIR \n");
			logDFTThreeModeSmoothExpand(EIRInit, origEIR, daysInYear, EIRX);
		}
		if(ifUseFC){
			printf("Calculating inverse discrete Fourier transform over full year \n");
			calcInverseDFTExp(EIRInit, daysInYear, FCEIR, FCEIRX);
		}
		if(ifrotateEIR){
			rotateArray(EIRInit, daysInYear, EIRRotateAngle);
			printf("Rotating EIR \n");
		}
	}
	else{
		convertLengthToFullYear(EIRInit, EIR);
	}
	if(ifprintEIR){
		PrintArray(fnametestentopar, origEIRname, origEIR, intervalsPerYear);
		PrintArray(fnametestentopar, shortEIRname, EIR, intervalsPerYear);
		PrintArray(fnametestentopar, longEIRname, EIRInit, daysInYear);
	}

	convertLengthToFullYear(humanInfectivityInit, initialKappa);
	if(ifprintKvi){
		PrintArray(fnametestentopar, shortKviname, initialKappa, intervalsPerYear);
		PrintArray(fnametestentopar, longKviname, humanInfectivityInit, daysInYear);
	}

  
    printvectors=0;
    if ( printvectors) {
        fWriteString(0, "Inteval = ");
        fWriteInt(0, interval);
        fWriteString(0, "\n");
        //end of fWrite
        fWriteString(0, "Intervals Per Year = ");
        fWriteInt(0, intervalsPerYear);
        fWriteString(0, "\n");
        //end of fWrite
        fWriteString(0, "Days in Year = ");
        fWriteInt(0, daysInYear);
        fWriteString(0, "\n");
        //end of fWrite
        for ( i=1;i<=intervalsPerYear; i++) {
            fWriteString(0, "initialKappa(");
            fWriteInt(0, i);
            fWriteString(0, ") = ");
            fWriteInt(0, initialKappa[i - 1]);
            fWriteString(0, "\n");
            //end of fWrite
        }
        for ( i=1;i<=daysInYear; i++) {
            fWriteString(0, "humanInfectivityInit(");
            fWriteInt(0, i);
            fWriteString(0, ")= ");
            fWriteInt(0, humanInfectivityInit[i - 1]);
            fWriteString(0, "\n");
            //end of fWrite
        }
        for ( i=1;i<=daysInYear; i++) {
            fWriteString(0, "EIR(");
            fWriteInt(0, i);
            fWriteString(0, ")= ");
            fWriteInt(0, EIR[i - 1]);
            fWriteString(0, "\n");
            //end of fWrite
        }
    }
    /*
     Set values of the initial estimate for the mosquito emergence
     rate. 
	 If we have already calcuated the mosquito emergence rate for the
	 given parameters separately, we can simply use that (and later
	 test the resulting EIR if it matches.
	 Otherwise, we just use a multiple of the EIR. The value of this 
	 vector may not be very important, but it may speed up the root
	 finding algorithm (but probably not).
	 (2008.10.20: It appears to make no difference to the speed.)
    */
	if(ifUseNv0Guess){
		// Read from file.
	}
	else{
		for ( i=1;i<=daysInYear; i++) {
			mosqEmergeRateInitEstimate[i - 1] = 
					EIRInit[i - 1]*popSizeInit*popSizeInit*hostAvailabilityRateInit;
		}
    }
    // Set all values of mosqEmergeRate to zero before calling C.
    for ( i=1;i<=daysInYear; i++) {
        mosqEmergeRate[i - 1]=0;
    }

        /*

     CALL C FUNCTION HERE !!!!!!!!!!!!!!!!!!!!!!!!!

    */
    CalcInitMosqEmergeRate(mosqEmergeRate, din, mosqRestDuration, 
		 EIPDuration, hti, mhti, popSizeInit, hostAvailabilityRateInit, 
		 mosqSeekingDeathRate, mosqSeekingDuration, mosqProbBiting, 
		 mosqProbFindRestSite, mosqProbResting, mosqProbOvipositing, 
		 humanInfectivityInit, EIRInit, mosqEmergeRateInitEstimate,
		 fnametestentopar);
}


void convertLengthToFullYear (double FullArray[daysInYear], double* ShortArray) {
    /*
    ***************************************************************!
    ******************* convertLengthToFullYear *******************!
    ***************************************************************!
     Note that ShortArray is assumed to be a pointer to a double array
	 of length intervalsPerYear. We do not explicitly check this.

     This subroutine converts vectors of length intervalsPerYear
     to daysInYear.
    
     For now, we assume that we will use intervalsPerYear and
     daysInYear as they are defined in global.f. We do not make this
     subroutine a general subroutine that converts from a given 
     length to another given length.
    
     Note that ShortArray should be of type real while FullArray
     should be of type real*8.
    
     Note: We also assume that:
     daysInYear = interval*intervalsPerYear.
    */

    int i;
    int j;
    /*
     TODO We assume that daysInYear = interval*intervalsPerYear
     Throw an error if this is not true.
    */
    for ( i=1;i<=intervalsPerYear; i++) {
        for ( j=1;j<=interval; j++) {
            FullArray[(i-1)*interval+j - 1]=(ShortArray[i - 1]);
        }
    }
}

void surfaceAreaAgeConversion () {

    double avbites_6;
    //agemin0 is the lower bound of the youngest age group used in this calculation
    double agemin0;
    int i;
    for ( i=1;i<=nages; i++) {
        btprop[i - 1]=bsa_prop[i - 1];
    }
    i=1;
    avbites_6=0.0;
    while( agemin[i - 1] <  6) {
        agemin0=agemin[i - 1];
        if ( agemin[i - 1] <  0.5 &&  agemax[i - 1] >  0.5) {
            agemin0=0.5;
        }
        if ( agemax[i - 1] >  0.5) {
            avbites_6=avbites_6+btprop[i - 1]*(agemax[i - 1]-agemin0)/5.5;
        }
        i=i+1;
    }
    biteratio_6=avbites_6/(1-avbites_6);
}

void inputEIR () {
    /*
    
    Reads in the estimates of the EIR for each village and each day
    and converts this into EIR estimates per five day period
    assuming that the annual cycle repeated during the pre-intervention period
    
    */


	// TODO: Properly include new simulationMode for entomological model. 
	// Depending on flags, this should either take in given EIR and smooth through
	// a DFT or directly take in Fourier coefficients and create an EIR over time.
    int mpcday;
    int j;
    double EIRdaily;
    //The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR;


	// We define these to be able to print out the original EIR.
	int ifprintorigEIR = 0;
	char origeirname[15] = "originalEIR";	

	int ifprintEIRaIDFT = 0;
	char eiraidftname[15] = "EIRafterIDFT";

    //initialise all the EIR arrays to 0
    if ( simulationMode !=  transientEIRknown) {
        for ( j=1;j<=intervalsPerYear; j++) {
            EIR[j - 1]=0.0;
            no[j - 1]=0.0;
        }
    }
    else {
        for ( j=1;j<=maxIntervals; j++) {
            intEIR[j - 1]=0.0;
            ino[j - 1]=0.0;
        }
    }
    minEIR=min_EIR_mult*averageEIR();
    mpcday=0;
    EIRdaily=get_eir_daily(mpcday);
    while( EIRdaily !=  missing_value) {
        updateEIR(mpcday, max(EIRdaily, minEIR));
        mpcday=mpcday+1;
        EIRdaily=get_eir_daily(mpcday);
    }
    //	Calculate total annual EIR 
    annualEIR=0.0;
    if ( simulationMode !=  transientEIRknown) {
        for ( j=1;j<=intervalsPerYear; j++) {
            annualEIR=annualEIR+interval*EIR[j - 1];
        }
    }
    else {
        annualEIR=-9.99;
    }
	
	// We copy EIR into origEIR
	for(j=0; j<EIRX; j++){
		origEIR[j] = EIR[j];
	}

	// For now we assume that we can manipulate EIR depending on the value of FTSmoothEIR.
	// We are assuming that the simulation mode is not set to transientEIRknown.
	// This needs to be rewritten properly once we have introduced a new simulationMode.
	// If FTSmoothEIR is 0, we do nothing.
	// If FTSmoothEIR is 1, we smooth the EIR using the first 3 modes of the discrete
	//		Fourier Transform
	if(FTSmoothEIR==1){
		// Print original EIR
		if(ifprintorigEIR){
			PrintArray(fnametestentopar, origeirname, origEIR, EIRX);
		}
		
		logDFTThreeModeSmooth(EIR, origEIR, EIRX);
	}
	if(ifUseFC==1){
		calcInverseDFTExp(EIR, EIRX, FCEIR, FCEIRX);
		if(ifprintEIRaIDFT == 1){
			PrintArray(fnametestentopar, eiraidftname, EIR, EIRX);
		}
	}


	if(ifrotateEIR){
		rotateArray(EIR, EIRX, EIRRotateAngle);
	}
}

void updateEIR (int day, double EIRdaily) {
    /*
    Processes each daily EIR estimate, allocating each day in turn to the appropriate 
    time period. 
    EIRdaily is the value of the daily EIR read in from the .XML file :
    
    */

    int i1;
    int istep;
    /*
    istep is the time period to which the day is assigned.  The result of the division is automatically
    rounded down to the next integer
    */
    istep=1+(day-1)/interval;
    if ( simulationMode !=  transientEIRknown) {
        i1=modIntervalsPerYear(istep);
        no[i1 - 1]=no[i1 - 1]+1.0;
        //EIR() is the arithmetic mean of the EIRs assigned to the 73 different recurring time points
        EIR[i1 - 1]=((EIR[i1 - 1]*(no[i1 - 1]-1))+EIRdaily)/no[i1 - 1];
    }
    else {
        i1=istep;
        ino[i1 - 1]=ino[i1 - 1]+1.0;
        intEIR[i1 - 1]=((intEIR[i1 - 1]*(ino[i1 - 1]-1))+EIRdaily)/ino[i1 - 1];
    }
}

void clearEntomologyParameters () {
    //Deallocate memory, clean up

    free(no);
    free(kappa);
    free(initialKappa);
    free(EIR);
    free(ino);
    free(intEIR);
	free(origEIR);
	free(FCEIR);
}

double averageEIR () {
    /*
    
    Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
    
    */

    int i;
    double valaverageEIR;
    i=0;
    valaverageEIR=0.0;
    while( get_eir_daily(i) !=  missing_value) {
        valaverageEIR=valaverageEIR+get_eir_daily(i);
        i=i+1;
    }
    valaverageEIR=valaverageEIR/i;
    return valaverageEIR;
}

double getBiteRatio (double ageyrs) {

    int i;
    double valgetBiteRatio;
    i=1;
    //60 yrs is the last cutpoint in the human growth curve
    while( agemax[i - 1] <  ageyrs &&  agemax[i - 1] <  60) {
        i=i+1;
    }
    /*
    the btprop() vector contains proportions of the total bites received by a host of this age 
    when competing with an adult 
    */
    valgetBiteRatio=btprop[i - 1]/(1-btprop[i - 1]);
    return valgetBiteRatio;
}

float hcalculate (float *cumEIR, float efficacy, float age_adj_EIR, float baseAvailToMos) {
    /*
    
    1.  Calculates h from the EIR measured on adults where
        h is the expected number of epidemiological inoculations
    2   Calculates the updated values of the pre-erythrocytic exposure and
        passes this back to the calling routine
    Requires the five-day EIR, adjusted for age as input. 
    cumEIR: is the pre-erythrocytic exposure;
    efficacy: Efficacy of a pre-erythrocytic vaccine
    
    */

    double normp;
    float S;
    float  ExpectedInfectionRate;
    float Infectionrate;
    //The age-adjusted EIR, possibly adjusted for bed nets.
    float effectiveEIR;
    float valhcalculate;
    if ( ITN) {
        effectiveEIR=age_adj_EIR*sqrt(Pu1/Pu0);
    }
    else {
        effectiveEIR=age_adj_EIR;
    }
    ExpectedInfectionRate=effectiveEIR*baseAvailToMos*susceptibility*interval;
    if ( isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
        Infectionrate=(W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
        valhcalculate=Infectionrate;
    }
    else if( isOptionIncluded(modelVersion, lognormalMassAction)) {
        normp=W_UNIFORM();
        Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
        valhcalculate=Infectionrate;
        //Bad (duplicated) code
    }
    else if( isOptionIncluded(modelVersion, lognormalMassActionPlusPreImm)) {
        normp=W_UNIFORM();
        Infectionrate=sampleFromLogNormal(normp, log(ExpectedInfectionRate)-0.5*pow(InfectionrateShapeParam, 2), InfectionrateShapeParam);
        S=(1.0+pow((*cumEIR/Xstar_p), gamma_p));
        S=Simm+(1.0-Simm)/S;
        S=S*(Sinf+(1-Sinf)/(1+effectiveEIR/Estar));
        valhcalculate=S*Infectionrate;
    }
    else {
        S=(1.0+pow((*cumEIR/Xstar_p), gamma_p));
        S=Simm+(1.0-Simm)/S;
        S=S*(Sinf+(1-Sinf)/(1+effectiveEIR/Estar));
        valhcalculate=S*effectiveEIR*interval;
        /*
        	  hcalculate=S*effectiveEIR*interval*baseAvailToMos 
               for heterogeneity
        */
    }
    //Introduce the effect of vaccination. Note that this does not affect cumEIR.
    if ( isOptionIncluded(vaccineType, preerythrocytic_reduces_h)) {
        valhcalculate=valhcalculate*(1-efficacy);
    }
    //Update pre-erythrocytic immunity
    *cumEIR =*cumEIR+interval*effectiveEIR;
    /*
          for heterogeneity	
          cumEIR=cumEIR+(interval*effectiveEIR*baseAvailToMos)
    */
    return valhcalculate;
}

double calculateEIR () {
    /*
    
    Calculates EIR (in adults), based on vectorial capacity or looks up EIR
    in the input data.
    
    time: Time since start of simulation .
    
    */

    double EIRadult;
    double Pc0;
    double Puz;
    double Pcz;
    double dz;
    double s0t;
    double szt;
    //	where the full model, with estimates of human mosquito transmission is in use, use this:
    double valcalculateEIR;
    switch (simulationMode) {
        case equilibriumMode:
        EIRadult=EIR[modIntervalsPerYear(t) - 1];
        break;
        case transientEIRknown:
        /*
        where the EIR for the intervention phase is known, obtain this from the
        intEIR array
        */
        EIRadult=intEIR[tstep - 1];
        break;
        case dynamicEIR:
        if ( tstep ==  1) {
            EIRadult=EIR[modIntervalsPerYear(t) - 1];
        }
        else {
            if ( ITN) {
                Pc0=pow(Pu0, c);
                Puz=Pu0-z*(Pu0-Pu1);
                Pcz=pow(Puz, c);
                dz=(1-Pu0)/(1-Puz);
                s0t=initialKappa[modIntervalsPerYear(t-nearbyint(nspore)) - 1]*Pc0/(1-Pu0);
                szt=kappa[modIntervalsPerYear(t-nearbyint(nspore)) - 1]*Pcz/(1-Puz);
                EIRadult=EIR[modIntervalsPerYear(t) - 1]*dz*szt/s0t;
            }
            else {
                EIRadult=EIR[modIntervalsPerYear(t) - 1]*kappa[modIntervalsPerYear(t-nearbyint(nspore)) - 1]/initialKappa[modIntervalsPerYear(t-nearbyint(nspore)) - 1];
            }
        }
        break;
    }
    valcalculateEIR=EIRadult;
    return valcalculateEIR;
}


void logDFTThreeModeSmooth(double* smoothArray, double* originalArray, int aLength) {
    /*
    ***************************************************************
    ******************** logDFTThreeModeSmooth ********************
    ***************************************************************
    *  Given a positive array, originalArray, of length aLength,
	*  this routine exponentiates the inverse discrete Fourier 
	* tranform of the first three modes of the natural logarithm of 
	* the array to smooth out the array to produce smoothArray.
    *
	* All elements of originalArray are assumed to be strictly
	* positive.
    *
	* smoothArray is an OUT parameter.
	* originalArray and aLength are IN parameters.
    */

    int t;
	// Period
	double P;
	// Frequency
	double w;
    
	// Fourier Coefficients
	double a0;
	double a1;
	double b1;
	double a2;
	double b2;

	double tempsuma0;
	double tempsuma1;
	double tempsumb1;
	double tempsuma2;
	double tempsumb2;

	double yt;		// Temporary log of array.

	// We define these to be able to print out smoothArray.
	int ifprintsmootharray = 0;
	char saname[15] = "SmoothArray";

	P=(double)aLength;
	w = 2*M_PI/P;

	// Calculate first three Fourier modes
    tempsuma0 = 0;
	tempsuma1 = 0;
	tempsumb1 = 0;
	tempsuma2 = 0;
	tempsumb2 = 0;

	for (t=1;t<=aLength; t++) {
		yt = log(originalArray[t-1]);
		tempsuma0 = tempsuma0+yt;
		tempsuma1 = tempsuma1 + (yt*cos(w*t));
		tempsumb1 = tempsumb1 + (yt*sin(w*t));
		tempsuma2 = tempsuma2 + (yt*cos(2*w*t));
		tempsumb2 = tempsumb2 + (yt*sin(2*w*t));       
    }
	a0 = (1/P)*tempsuma0;
	a1 = (2/P)*tempsuma1;
	b1 = (2/P)*tempsumb1;
	a2 = (2/P)*tempsuma2;
	b2 = (2/P)*tempsumb2;

	// Calculate inverse discrete Fourier transform
	for (t=1; t<=aLength; t++){
		smoothArray[t-1] = 
			exp(a0 + a1*cos(w*t) + b1*sin(w*t) + a2*cos(2*w*t) + b2*sin(2*w*t));
	}

	if(ifprintsmootharray){
			PrintArray(fnametestentopar, saname, smoothArray, aLength);

			printf("a0 = %f; \n", a0);
			printf("a1 = %f; \n", a1);
			printf("b1 = %f; \n", b1);
			printf("a2 = %f; \n", a2);
			printf("b2 = %f; \n", b2);
			// getchar();

		}


}



void logDFTThreeModeSmoothExpand(double* smoothArray, 
								 double* originalArray, int SALength, int OALength) {
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

    int t;
	// Period
	double Poa;
	double Psa;
	// Frequency
	double woa;
	double wsa;
    
	// Fourier Coefficients
	double a0;
	double a1;
	double b1;
	double a2;
	double b2;

	double tempsuma0;
	double tempsuma1;
	double tempsumb1;
	double tempsuma2;
	double tempsumb2;

	double yt;		// Temporary log of array.

	// We define these to be able to print out smoothArray.
	int ifprintsmootharray = 0;
	char saname[15] = "SmoothArray";

	Poa=(double)OALength;
	Psa=(double)SALength;

	woa = 2*M_PI/Poa;
	wsa = 2*M_PI/Psa;

	// Calculate first three Fourier modes
    tempsuma0 = 0;
	tempsuma1 = 0;
	tempsumb1 = 0;
	tempsuma2 = 0;
	tempsumb2 = 0;

	for (t=1;t<=OALength; t++) {
		yt = log(originalArray[t-1]);
		tempsuma0 = tempsuma0+yt;
		tempsuma1 = tempsuma1 + (yt*cos(woa*t));
		tempsumb1 = tempsumb1 + (yt*sin(woa*t));
		tempsuma2 = tempsuma2 + (yt*cos(2*woa*t));
		tempsumb2 = tempsumb2 + (yt*sin(2*woa*t));       
    }
	a0 = (1/Poa)*tempsuma0;
	a1 = (2/Poa)*tempsuma1;
	b1 = (2/Poa)*tempsumb1;
	a2 = (2/Poa)*tempsuma2;
	b2 = (2/Poa)*tempsumb2;

	// Calculate inverse discrete Fourier transform
	for (t=1; t<=SALength; t++){
		smoothArray[t-1] = 
			exp(a0 + a1*cos(wsa*t) + b1*sin(wsa*t) + 
			    a2*cos(2*wsa*t) + b2*sin(2*wsa*t));
	}

	if(ifprintsmootharray){
			PrintArray(fnametestentopar, saname, smoothArray, SALength);
		}


}



void calcInverseDFTExp(double* tArray, int aL, double* FC, int FCL) {
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

    int t;
	int n;
	// Period
	double P;
	// Frequency
	double w;
	// Number of Fourier Modes.
	int Fn;
	
	double temp;
    
	P=(double)aL;
	w = 2*M_PI/P;


	if((FCL%2)==0){
		printf("The number of Fourier coefficents should be odd.\n");
		getchar();
	}

	Fn = (FCL-1)/2;

	// Calculate inverse discrete Fourier transform
	for (t=1; t<=aL; t++){
		temp = FC[0];
		if(Fn>0){
			for(n=1;n<=Fn;n++){
				temp = temp + FC[2*n-1]*cos(n*w*t) + FC[2*n]*sin(n*w*t);
			}
		}
		tArray[t-1] = exp(temp);
	}
}

void rotateArray(double* rArray, int aLength, double rAngle) {
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

    int i;
	int rotIndex;
	double* tempArray;
	int tIndex;


	// We define these to be able to print out smoothArray.
	int ifprintrarray = 1;
	char raprename[20] = "PrerotationArray";
	char rapostname[20] = "PostrotationArray";

	tempArray = (double *)malloc((aLength)*sizeof(double));

	rotIndex = (int)((rAngle*aLength)/(2*M_PI));


	if(ifprintrarray){
		PrintArray(fnametestentopar, raprename, rArray, aLength);
	}


	for (i=0;i<aLength; i++) {
		tIndex = (i+rotIndex) % aLength;
		tempArray[tIndex] = rArray[i];
    }
	

	for (i=0; i<aLength; i++){
		rArray[i] = tempArray[i];
	}

	if(ifprintrarray){
		PrintArray(fnametestentopar, rapostname, rArray, aLength);
	}
	free(tempArray);
}

