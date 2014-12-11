/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#include "WithinHost/Infection/DescriptiveInfection.h"
#include "util/random.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <sstream>
#include <string>
#include <cmath>
#include <fstream>

namespace OM {
namespace WithinHost {
using namespace util;

// static class variables (see description in header file):
double DescriptiveInfection::meanLogParasiteCount[numDurations][numDurations];
double DescriptiveInfection::sigma0sq;
double DescriptiveInfection::xNuStar;

bool bugfix_max_dens = true, bugfix_innate_max_dens = true;


// ———  static init/clear ———

void DescriptiveInfection::init (const Parameters& parameters) {
    // Error checks
    if( sim::oneTS().inDays() != 5 ){
        // To support non-5-day time-step models, either different data would
        // be needed or times need to be adjusted when accessing
        // meanLogParasiteCount. Probably the rest would be fine.
        throw util::xml_scenario_error ("DescriptiveInfection only supports using an interval of 5");
    }
    // Bug fixes: these are enabled by default but may be off in old parameterisations
    bugfix_innate_max_dens = util::ModelOptions::option (util::INNATE_MAX_DENS);
    // Warning: if MAX_DENS_CORRECTION is off, infections not yet at the blood
    // stage could result in BSVEfficacy and potentially innateImmSurvFact
    // being applied to timeStepMaxDensity more than once in some cases.
    bugfix_max_dens = util::ModelOptions::option (util::MAX_DENS_CORRECTION);
    
    // Read parameters
    sigma0sq=parameters[Parameters::SIGMA0_SQ];
    xNuStar=parameters[Parameters::X_NU_STAR];
    
    // Read file empirical parasite densities
    string densities_filename = util::CommandLine::lookupResource ("densities.csv");
    ifstream f_MTherapyDensities( densities_filename.c_str() );
    if( !f_MTherapyDensities.good() ){
        throw util::base_exception( string("Cannot read ").append(densities_filename), util::Error::FileIO );
    }
    
    //read header of file (unused)
    string csvLine;
    getline(f_MTherapyDensities,csvLine);

    // read every line from the stream
    while (getline(f_MTherapyDensities, csvLine)) {
        //Empirical description of single Malaria infections in naive individuals
        //counter variables, i stands for 5 day time interval, j for duration of infection
        int i;
        int j;
        double meanlogdens;

        std::istringstream csvStream(csvLine);

        string csvField1, csvField2, csvField3;

        // read every element from the line that is seperated by commas
        getline(csvStream, csvField1, ',');
        getline(csvStream, csvField2, ',');
        getline(csvStream, csvField3, ',');

        istringstream csvNum1(csvField1), csvNum2(csvField2), csvNum3(csvField3);

        csvNum1 >> i;
        csvNum2 >> j;
        csvNum3 >> meanlogdens;
        
        if( i < 1 || i > numDurations || j < 1 || j > numDurations ){
            throw util::base_exception( string("error in data file: ")
                .append(densities_filename), Error::InputResource );
        }

        //fill initial matrix
        meanLogParasiteCount[i-1][j-1]=meanlogdens;
        //fill also the triangle that will not be used (to ensure everything is initialised)
        if (j!=i) {
            meanLogParasiteCount[j-1][i-1]=0.0;
        }

    }
}


// ———  non-static init/destruction  ———

DescriptiveInfection::DescriptiveInfection () :
        m_duration(infectionDuration()),
        notPrintedMDWarning(true)
{
    assert( sim::oneTS().inDays() == 5 );
}

SimTime DescriptiveInfection::infectionDuration() {
    // Forgive the excess precision; it just avoids having to update all expected results
    double dur_mean = 5.1300001144409179688;
    double dur_sigma = 0.80000001192092895508;
    double dur=random::log_normal(dur_mean, dur_sigma);
    
    //TODO:
    // Model did say infection is cleared on day dur+1 converted to a time-step
    // let interval = sim::oneTS().inDays() in:
    // ((1+floor(dur))/interval); now it says the last interval is:
    // floor((1+dur)/interval)-1 = floor((dur+1-interval)/interval)
    // Is this reasonable, or should we change?
    return sim::fromDays( static_cast<int>(1.0 + dur) ) - sim::oneTS();
}


// ———  time-step updates  ———

void DescriptiveInfection::determineDensities(double ageInYears,
                                              double cumulativeh,
                                              double cumulativeY,
                                              double &timeStepMaxDensity,
                                              double innateImmSurvFact,
                                              double bsvFactor)
{
    // Age of patent blood stage infection. Note: liver stage is fixed at one
    // 5-day time step and prepatent blood stage is latentp - 1 time steps.
    SimTime infage = sim::ts0() - m_startDate - latentP;
    if ( infage < sim::zero()) {
        m_density = 0.0;
        if (bugfix_max_dens) timeStepMaxDensity = 0.0;
    }else{
        timeStepMaxDensity = 0.0;
        
        int32_t infAge = min( infage.inSteps(), maxDurationTS );
        int32_t infDur = min( m_duration.inSteps(), maxDurationTS );
        m_density=max (exp(meanLogParasiteCount[infAge][infDur]), 1.0);
        
        // The expected parasite density in the non naive host (AJTM p.9 eq. 9)
        // Note that in published and current implementations Dx is zero.
        m_density = pow(m_density, immunitySurvivalFactor(ageInYears, cumulativeh, cumulativeY));
        
        //Perturb m_density using a lognormal
        double varlog = sigma0sq / (1.0 + (cumulativeh / xNuStar));
        double stdlog = sqrt(varlog);
        
        /*
        This code samples from a log normal distribution with mean equal to the predicted density
        n.b. AJTM p.9 eq 9 implies that we sample the log of the density from a normal with mean equal to
        the log of the predicted density.  If we really did the latter then this bias correction is not needed.
        */
        double meanlog = log(m_density) - stdlog*stdlog / 2.0;
        if (stdlog > 0.0000001) {
            // Calculate the expected density on the day of sampling:
            m_density = random::log_normal(meanlog, stdlog);
            // Calculate additional samples for T-1 days (T=sim::oneTS().inDays()):
            if( true /*was sim::oneTS().inDays() > 1, which is always true in this model*/ ){
                double normp = pow( random::uniform_01(), 1.0 / (sim::oneTS().inDays() - 1) );
                /*
                To mimic sampling T-1 repeated values, we transform the sampling
                distribution and use only one sampled value, which has the sampling
                distribution of the maximum of T-1 values sampled from a uniform.
                The probability density function of this sampled random var is distributed
                according to a skewed distribution (defined in [0,1]) where the
                exponent (1/(T-1)) arises because each of T-1 sampled
                values would have this probability of being the maximum.
                */
                timeStepMaxDensity = random::sampleFromLogNormal(normp, meanlog, stdlog);
            }
            timeStepMaxDensity = std::max(m_density, timeStepMaxDensity);
        }
        if (timeStepMaxDensity > maxDens && notPrintedMDWarning){
            cerr << "TSMD hit limit:\t" << m_density << ",\t" << timeStepMaxDensity << endl;
            notPrintedMDWarning = false;
        }
        m_density = min(m_density, maxDens);
        timeStepMaxDensity = min(timeStepMaxDensity, maxDens);
        
        //Compute the proportion of parasites remaining after innate blood stage effect
        m_density *= innateImmSurvFact;
        
        //Include here the effect of blood stage vaccination
        m_density *= bsvFactor;
        
        m_cumulativeExposureJ += sim::oneTS().inDays() * m_density;
    }
    
    if (bugfix_innate_max_dens) timeStepMaxDensity *= innateImmSurvFact;
    timeStepMaxDensity *= bsvFactor;
}


// ———  checkpointing  ———

DescriptiveInfection::DescriptiveInfection (istream& stream) :
        Infection(stream), m_duration(sim::never())
{
    m_duration & stream;
    notPrintedMDWarning & stream;
}
void DescriptiveInfection::checkpoint (ostream& stream) {
    Infection::checkpoint (stream);
    m_duration & stream;
    notPrintedMDWarning & stream;
}

}
}
