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

int MolineauxInfection::delta;

double MolineauxInfection::C;
double MolineauxInfection::sigma;
double MolineauxInfection::rho;
double MolineauxInfection::beta;
double MolineauxInfection::sProb;
double MolineauxInfection::q;
double MolineauxInfection::mu_m;
double MolineauxInfection::sigma_m;
double MolineauxInfection::k_c;
double MolineauxInfection::k_m;
double MolineauxInfection::Pstar_v;
double MolineauxInfection::kappa_c;
double MolineauxInfection::kappa_v;
double MolineauxInfection::kappa_m;

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

	C=1.0, sigma=0.02,rho=0, beta=0.01, sProb=0.02, q=0.3, mu_m=16.0, sigma_m=10.4,
	k_c=0.2, k_m=0.04, Pstar_v=30.0, kappa_c=3.0, kappa_v=3.0, kappa_m=1.0;

	delta = 8;
}

MolineauxInfection::MolineauxInfection(uint32_t protID):
		CommonInfection(protID)
{
	for(int i=0;i<v; i++)
	{
		P[i] = 0.0;
		m[i] = 0.0;
		growthRate[i] = 0.0;

		while(m[i]<1.0)
			m[i]=random::gauss(mu_m, sigma_m);

		for(int tau=0; tau<taus; tau++)
			laggedP[tau][i] =  0.0;

		variantSpecificSummation[i]=0.0;
	}
	for(int tau=0; tau<taus;tau++)
		laggedPc[tau] = 0.0;

	P[0] = 0.1;
	variantTranscendingSummation = 0.0;
	Pstar_c = k_c*pow(random::gauss(Params::MEAN_LOCAL_MAX_DENSITY, Params::SD_LOCAL_MAX_DENSITY),10.0);
	Pstar_m = k_m*pow(random::gauss(Params::MEAN_DIFF_POS_DAYS, Params::SD_DIFF_POS_DAYS), 10.0);
}

void MolineauxInfection::updateGrowthRateMultiplier(){

	double Sc = 1/(1 + pow(_density/Pstar_c, kappa_c));
	double Sm = ((1-beta)/(1+pow(getVariantTranscendingSummation()/Pstar_m, kappa_m)))+beta;
	double S[v];

	double sigma_Qi_Si=0.0;
	double sigma_S=0.0;
	double sigmaP = 0.0;

	for(int j=0;j<v; j++)
		sigmaP+=P[j];

	for(int i=0; i<v; i++)
	{
		S[i] = 1/(1+pow(getVariantSpecificSummation(i,P[i])/Pstar_v, kappa_v));
		sigma_Qi_Si+= pow(q, (double)i+1)*S[i];
		sigma_S+=S[i];
	}

	for(int i=0;i<v;i++)
	{
		double p_i;
		if(S[i]<0.1)
			p_i = 0.0;
		else
			p_i = pow(q, (double)i+1)*S[i]/sigma_Qi_Si;

		double newPi;
		newPi = ((1-sProb) * P[i]+sProb*p_i*sigmaP)*m[i]*S[i]*Sc*Sm;
		if(newPi<1.0e-5)
			newPi = 0.0;

		if(P[i]==0)
		{
			growthRate[i] = 1;
			P[i] = newPi;
		}
		else
			growthRate[i] = sqrt(newPi/P[i]);
	}
}

bool MolineauxInfection::updateDensity(double survivalFactor, int ageOfInfection){
	double newDensity = 0.0;

	if(ageOfInfection == 0)
		_density = P[0];
	else
	{
		for(int i=0;i<v;i++)
		{
			P[i] *= growthRate[i] * survivalFactor;

			if(P[i]<1.0e-5)
				P[i] = 0.0;
			else
				newDensity += P[i];
		}

		_density = newDensity;
	}

	_cumulativeExposureJ += Global::interval * _density;

	if(_density>1.0e-5)
	{
		if(ageOfInfection%2==0)
		{
			updateGrowthRateMultiplier();
		}
		return false;
	}
	else
		return true;
}

double MolineauxInfection::getVariantSpecificSummation(int i, double Pcurrent){
//The effective exposure is computed by adding in the 8-day lagged parasite density (i.e. 4 time steps)
//and decaying the previous value for the effective exposure with decay parameter 2*sigma (the 2 arises because
//the time steps are two days and the dimension of sigma is per day.

	variantSpecificSummation[i] = (variantSpecificSummation[i] * exp(-2.0*sigma))+laggedP[3][i];
	for(int tau=1; tau<taus; tau++)
		laggedP[4-tau][i] = laggedP[3-tau][i];

	laggedP[0][i] = Pcurrent;

	return variantSpecificSummation[i];
}

double MolineauxInfection::getVariantTranscendingSummation(){

	variantTranscendingSummation = (variantTranscendingSummation * exp(-2.0*rho))+laggedPc[3];
	for (int tau=1; tau<taus; tau++)
		laggedPc[4-tau] = laggedPc[3-tau];

	if(_density < C)
		laggedPc[0] = _density;
	else
		laggedPc[0] =C;

	return variantTranscendingSummation;
}

MolineauxInfection::MolineauxInfection (istream& stream) :
    CommonInfection(stream)
{
    variantTranscendingSummation & stream;
    for(int i=0;i<v;i++){
	m[i] & stream;
    	growthRate[i] & stream;
	P[i] & stream;
	variantSpecificSummation[i] & stream;
	for(int j=0;j<taus;j++)
	    laggedP[j][i] & stream;
    }
    for(int j=0;j<taus;j++)
	laggedPc[j] & stream;
    Pstar_c & stream;
    Pstar_m & stream;
}

void MolineauxInfection::checkpoint (ostream& stream) {
    CommonInfection::checkpoint (stream);

    variantTranscendingSummation & stream;
    for(int i=0;i<v;i++){
	m[i] & stream;
    	growthRate[i] & stream;
	P[i] & stream;
	variantSpecificSummation[i] & stream;
	for(int j=0;j<taus;j++)
	    laggedP[j][i] & stream;
    }
    for(int j=0;j<taus;j++)
	laggedPc[j] & stream;
    Pstar_c & stream;
    Pstar_m & stream;
}

}
}
