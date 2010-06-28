/*
This file is part of OpenMalaria.

Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/CommonWithinHost.h"
#include "inputData.h"
#include "util/random.h"
#include "util/errors.hpp"
#include "util/CommandLine.hpp"
#include "util/ModelOptions.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>



namespace OM {
namespace WithinHost {

using namespace OM::util;

double MolineauxInfection::mean_first_local_max;
double MolineauxInfection::sd_first_local_max;
double MolineauxInfection::mean_diff_pos_days;
double MolineauxInfection::sd_diff_pos_days;

CommonInfection* createMolineauxInfection (uint32_t protID) {
    return new MolineauxInfection (protID);
}

CommonInfection* checkpointedMolineauxInfection (istream& stream) {
    return new MolineauxInfection (stream);
}

void MolineauxInfection::initParameters(){
	if (Global::interval != 1)
		throw util::xml_scenario_error ("MolineauxInfection only supports scenarii using an interval of 1");

	CommonWithinHost::createInfection = &createMolineauxInfection;
	CommonWithinHost::checkpointedInfection = &checkpointedMolineauxInfection;

	mean_first_local_max = InputData.getParameter(Params::MEAN_LOCAL_MAX_DENSITY);
	sd_first_local_max = InputData.getParameter(Params::SD_LOCAL_MAX_DENSITY);
	mean_diff_pos_days = InputData.getParameter(Params::MEAN_DIFF_POS_DAYS);
	sd_diff_pos_days = InputData.getParameter(Params::SD_DIFF_POS_DAYS);
}

MolineauxInfection::MolineauxInfection(uint32_t protID):
		CommonInfection(protID)
{
	for(int i=0;i<v; i++)
	{
		P[i] = 0.0;
		m[i] = 0.0;
		growthRate[i] = 0.0;
		initP[i] = 0.0;

		// Molineaux paper, equation 11
		while(m[i]<1.0)
		{
		    m[i]=random::gauss(mu_m, sigma_m);
		}

		for(int tau=0; tau<taus; tau++)
		{
		    laggedP[tau][i] =  0.0;
		}

		variantSpecificSummation[i]=0.0;
	}
	for(int tau=0; tau<taus;tau++)
	{
		laggedPc[tau] = 0.0;
	}

	// the initial density is set to 0.1... The first chosen variant is the variant 1
	/*initP[0] =*/ P[0] = 0.1;
	variantTranscendingSummation = 0.0;

	// TODO: use static variables
	Pstar_c = k_c*pow(random::gauss(mean_first_local_max,sd_first_local_max),10.0);
	Pstar_m = k_m*pow(random::gauss(mean_diff_pos_days,sd_diff_pos_days),10.0);
}

void MolineauxInfection::updateGrowthRateMultiplier(){
        // The immune responses are represented by the variables
        // Sc (probability that a parasite escapes control by innate and variant-transcending immune response)
        // Sm (                        "                      acquired and variant-transcending immune response)
        // S[i] (                      "                      acquired and variant-specific immune response)
	double Sc = 1.0/(1.0 + pow(_density/Pstar_c, kappa_c));

	//double Sm = ((1-beta)/(1+pow(getVariantTranscendingSummation()/Pstar_m, kappa_m)))+beta
	//optimization: Since kappa_m = 1, we don't use pow.
	double Sm = ((1.0-beta)/(1.0+(getVariantTranscendingSummation()/Pstar_m)))+beta;
	double S[v];

	double sigma_Qi_Si=0.0;
	
	for(int i=0; i<v; i++)
	{
		S[i] = 1.0/(1.0 + pow(getVariantSpecificSummation(i,P[i]) / Pstar_v, kappa_v));
		sigma_Qi_Si+= pow(q, (double)(i+1))*S[i];
		initP[i] = 0.0;
	}

	for(int i=0;i<v;i++)
	{
		// Molineaux paper equation 4
		// p_i: variant selection probability
		double p_i;
		if( S[i]<0.1 ){
		    p_i = 0.0;
		} else {
		    p_i = pow(q, (double)(i+1))*S[i]/sigma_Qi_Si;
		}

		// Molineaux paper equation 1
		// newPi: Variant density at t = t + 2
		// new variant density  = (the amount of this variant's parasites
		// which will not switch to another variant + the ones from other
		// variants switching to this variant) * this variant multiplication
		// factor * the probability that the parasites escape control by immune response.
		double newPi = ( (1.0-sProb) * P[i] + sProb*p_i*_density )*m[i]*S[i]*Sc*Sm;

		// Molineaux paper equation 2
		if(newPi<1.0e-5)
		{
		    newPi = 0.0;
		}

		// if P[i] == 0 then that means this variant wasn't expressed yet
		// or is extinct. If this variant is emerging in (t+2) the new variant
		// density is stored in the initP array,  so that we are able to add
		// the survival factor's effect to the emerging variant's density.
		if(P[i]==0) {
		    initP[i] = newPi;
		    growthRate[i] = 0.0;
		}
		else {
		    // important for checkpoints: leave initP[i] non-zero
		    growthRate[i] = sqrt(newPi/P[i]);
		}
	}
}

bool MolineauxInfection::updateDensity(double survivalFactor, int ageOfInfection){
	if(ageOfInfection == 0)
	{
	    _density = P[0];
	}
	else
	{
	    double newDensity = 0.0;

	    for(int i=0;i<v;i++)
	    {
		    // growthRate:
		    // p(t+1) = p(t) * sqrt(p(t+2)/p(t))
		    // p(t+2) = p(t+1) * sqrt(p(t+2)/p(t))
		    // p(t+2) = p(t) * sqrt(p(t+2)/p(t))^2...
		    P[i] *= growthRate[i];

		    // survivalFactor: effects of drugs, immunity and vaccines
		    P[i] *= survivalFactor;
		    initP[i] *= survivalFactor;

		    // if t+2: The new variant is now expressed. For already extinct
		    // variants this doesn't matter, since initP[i] = 0 for those variants.
		    if(P[i]==0 && ageOfInfection%2==0)
		    {
			P[i] = initP[i];
		    }

		    // The variant is extinct when variant's density < 1.0e-5
		    if(P[i]<1.0e-5)
		    {
			P[i] = 0.0;
		    }
		    else
		    {
			// Molineaux paper equation 3
			newDensity += P[i];
		    }
	    }

	    _density = newDensity;
	}

	_cumulativeExposureJ += Global::interval * _density;

	if(_density>1.0e-5)
	{
	    // if the infection isn't extinct and t = t+2
	    // then the growthRateMultiplier is adapted for t+3 and t+4
	    if(ageOfInfection%2==0)
	    {
		updateGrowthRateMultiplier();
	    }
	    return false;
	}
	else
	{
	    return true;
	}
}

double MolineauxInfection::getVariantSpecificSummation(int i, double Pcurrent){
    //The effective exposure is computed by adding in the 8-day lagged parasite density (i.e. 4 time steps)
    //and decaying the previous value for the effective exposure with decay parameter 2*sigma (the 2 arises because
    //the time steps are two days and the dimension of sigma is per day.

    //Molineaux paper equation 6
    size_t index = (Global::simulationTime % 8)/2;	// 8 days ago has same index as today
    variantSpecificSummation[i] = (variantSpecificSummation[i] * exp(-2.0*sigma))+laggedP[index][i];
    laggedP[index][i] = Pcurrent;

    return variantSpecificSummation[i];
}

double MolineauxInfection::getVariantTranscendingSummation(){

    //Molineaux paper equation 5
    size_t index = (Global::simulationTime % 8)/2;	// 8 days ago has same index as today
    variantTranscendingSummation = (variantTranscendingSummation * exp(-2.0*rho))+laggedPc[index];

    //Molineaux paper equation 8
    //We could use min here, but it seems that min has problems with static const double C
    laggedPc[index] = _density < C ? _density:C;

    return variantTranscendingSummation;
}

MolineauxInfection::MolineauxInfection (istream& stream) :
    CommonInfection(stream)
{
    variantTranscendingSummation & stream;
    for(int i=0;i<v;i++){
	m[i] & stream;
	initP[i] & stream;
	if( true /*initP[i] != 0.0*/ ){	// otherwise everything below is zero (checkpoint size optimization)
	    growthRate[i] & stream;
	    P[i] & stream;
	    variantSpecificSummation[i] & stream;
	    for(int j=0;j<taus;j++)
	    {
		laggedP[j][i] & stream;
	    }
	}
    }
    for(int j=0;j<taus;j++)
    {
	laggedPc[j] & stream;
    }
    Pstar_c & stream;
    Pstar_m & stream;
}

void MolineauxInfection::checkpoint (ostream& stream) {
    CommonInfection::checkpoint (stream);

    variantTranscendingSummation & stream;
    for(int i=0;i<v;i++){
    	m[i] & stream;
    	initP[i] & stream;
	if( true /*initP[i] != 0.0*/ ){
	    growthRate[i] & stream;
	    P[i] & stream;
	    variantSpecificSummation[i] & stream;
	    for(int j=0;j<taus;j++)
	    {
		laggedP[j][i] & stream;
	    }
	} else {
	    growthRate[i] = 0.0;
	    P[i] = 0.0;
	    variantSpecificSummation[i] = 0.0;
	    for(int j=0;j<taus;j++)
	    {
		laggedP[j][i] = 0.0;
	    }
	}
    }
    for(int j=0;j<taus;j++)
    {
	laggedPc[j] & stream;
    }
    Pstar_c & stream;
    Pstar_m & stream;
}

}
}
