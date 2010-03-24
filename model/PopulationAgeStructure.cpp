/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "PopulationAgeStructure.hpp"
#include "Global.h"
#include "inputData.h"
#include "util/errors.hpp"
#include "util/gsl.h"

#include <cmath>

namespace OM {
    using namespace OM::util;

#ifdef OM_HAVE_GSL
double AgeStructure::ageGroupBounds[ngroups+1];
double AgeStructure::ageGroupPercent[ngroups];

double AgeStructure::M1[ngroups];
double AgeStructure::M2[ngroups];
double AgeStructure::M[ngroups];
double AgeStructure::pred[ngroups];

double AgeStructure::mu0;
double AgeStructure::mu1;
double AgeStructure::alpha0;
double AgeStructure::alpha1;
double AgeStructure::rho;
#endif

int AgeStructure::maxTimestepsPerLife;
vector<double> AgeStructure::cumAgeProp;


void AgeStructure::init () {
    maxTimestepsPerLife = maxLifetimeDays / Global::interval;
    cumAgeProp.resize (maxTimestepsPerLife);
    
#ifdef OM_HAVE_GSL
    estimateRemovalRates();
    calcCumAgeProp();
#endif	// TODO: alternative load-from-XML functionality
}

int AgeStructure::targetCumPop (int ageTSteps, int targetPop)
{
    return (int) floor (cumAgeProp[maxTimestepsPerLife-1-ageTSteps] * targetPop + 0.5);
}


#ifdef OM_HAVE_GSL
void AgeStructure::estimateRemovalRates ()
{
    // mu1, alpha1: These are estimated here
    // mu0, alpha0: These are fixed alpha0=4.0, mu0 calculated from the other parameters
    // rho - population growth rate (input)
    
    //Get lower and upper age bounds for age groups and cumulative precentage of population from field data
    double sumperc = 0.0;
    const scnXml::AgeGroupPerC::GroupSequence& group = InputData().getDemography().getAgeGroup().getGroup();
    if (group.size() < ngroups - 1) {
	ostringstream msg;
	msg << "expected " << ngroups - 1 << " elements of \"group\" in demography->ageGroup (in scenario.xml)";
	throw util::xml_scenario_error (msg.str());
    }
    //Add age group for first month of life
    ageGroupBounds[0] = 0.0;
    ageGroupBounds[1] = 1.0 / 12.0;
    ageGroupPercent[0] = 0.0;
    for (int i = 1;i < ngroups; i++) {
	ageGroupBounds[i+1] = group[i-1].getUpperbound();
	ageGroupPercent[i] = group[i-1].getPoppercent();
	sumperc += ageGroupPercent[i];
    }
    sumperc = 100.0 / sumperc; // multiplier to get percentages
    for (int i = 0;i < ngroups; i++) {
	ageGroupPercent[i]  = ageGroupPercent[i] * sumperc;
    }
    /*
    RSS between observed and predicted log percentage of population in age groups
    is minimised for values of mu1 and alpha1
    calls setDemoParameters to calculate the RSS
    */
    /* NOTE: unused --- why?
    double tol = 0.00000000001;
    int maxcal = 500000;
    int npar = 2;
    int iw = 3; */
    double p1 = 0.371626412;
    double p2 = 0.841209593;
    // returns "double rss":
    gsl::minimizeCalc_rss (&setDemoParameters, p1, p2);
}

// Static method used by estimateRemovalRates
double AgeStructure::setDemoParameters (double param1, double param2)
{
    rho = 0.0;
    if (InputData().getDemography().getGrowthRate().present())
    	rho = InputData().getDemography().getGrowthRate().get();

    rho = rho * (0.01 * Global::yearsPerInterval);
    if (rho != 0.0)
	// Issue: in this case the total population size differs from populationSize,
	// however, some code currently uses this as the total population size.
	throw util::xml_scenario_error ("Population growth rate provided.");
    
    const double IMR = 0.1;
    double M_inf = -log (1 - IMR);
    
    mu1 = exp (param1) / 100;
    alpha1 = exp (param2) / 100;
    alpha0 = 4.0;
    mu0 = (M_inf - mu1 * (exp (alpha1 * 0.5) - 1) * alpha0) /
    (alpha1 * (1 - exp (-alpha0 * 0.5)));
    
    double sumpred = 0.0;
    for (int i = 0; i < ngroups - 1; i++) {
	double midpt = (ageGroupBounds[i+1] + ageGroupBounds[i]) * 0.5;
	M1[i] = mu0 * (1.0 - exp (-alpha0 * midpt)) / alpha0;
	M2[i] = mu1 * (exp (alpha1 * midpt) - 1.0) / alpha1;
	M[i]  = M1[i] + M2[i];
	pred[i] = (ageGroupBounds[i+1] - ageGroupBounds[i])
	* exp (-rho * midpt - M[i]);
	sumpred += pred[i];
    }
    for (int i = 0; i < ngroups - 1; i++) {
	pred[i] = pred[i] / sumpred * 100.0;
    }
    double L_inf = exp (-rho * 0.5 - M[1]);
    double M_nn  = -log (1.0 - 0.4 * (1 - exp (-M[1])));
    double L1    = 1.0 / 12.0 * exp (-rho / 24.0 - M_nn);
    double perc_inf = ageGroupPercent[0] + ageGroupPercent[1];
    ageGroupPercent[0] = perc_inf * L1 / L_inf;
    ageGroupPercent[1] = perc_inf - ageGroupPercent[0];
    
    double valsetDemoParameters = 0.0;
    for (int i = 0; i < ngroups - 1; i++) {
	double residual = log (pred[i]) - log (ageGroupPercent[i]);
	valsetDemoParameters += residual * residual;
    }
    return valsetDemoParameters;
}

void AgeStructure::calcCumAgeProp ()
{
    cumAgeProp[0] = 0.0;
    for (int j = 1;j < maxTimestepsPerLife; j++) {
	double ageYears = (maxTimestepsPerLife - j - 1) * Global::yearsPerInterval;
	double M1s = (mu0 * (1.0 - exp (-alpha0 * ageYears)) / alpha0);
	double M2s = (mu1 * (exp (alpha1 * ageYears) - 1.0) / alpha1);
	double Ms = M1s + M2s;
	double predperc = exp (-rho * ageYears - Ms);
	if (j < maxTimestepsPerLife - Global::maxAgeIntervals) {
	    predperc = 0.0;
	}
	cumAgeProp[j] = cumAgeProp[j-1] + predperc;
    }
    double totalCumPC = cumAgeProp[maxTimestepsPerLife-1];
    for (int j = 1;j < maxTimestepsPerLife; j++) {
	//Scale using the total cumAgeProp
	cumAgeProp[j] = cumAgeProp[j] / totalCumPC;
    }
}
#endif

}
