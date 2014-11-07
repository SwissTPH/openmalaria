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
#include "Transmission/TransmissionModel.h"
#include "Transmission/NonVectorModel.h"
#include "Transmission/VectorModel.h"
#include "Transmission/PerHost.h"

#include "Population.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Genotypes.h"
#include "Monitoring/Continuous.h"
#include "util/BoincWrapper.h"
#include "util/StreamValidator.h"
#include "util/CommandLine.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"

#include <cmath>
#include <cfloat>
#include <gsl/gsl_vector.h>

namespace OM { namespace Transmission {
namespace vectors = util::vectors;

TransmissionModel* TransmissionModel::createTransmissionModel (const scnXml::Entomology& entoData, int populationSize) {
  // Entomology contains either a list of at least one anopheles or a list of at
  // least one EIRDaily.
  const scnXml::Entomology::VectorOptional& vectorData = entoData.getVector();

  TransmissionModel *model;
  if (vectorData.present())
    model = new VectorModel(entoData, vectorData.get(), populationSize);
  else {
      const scnXml::Entomology::NonVectorOptional& nonVectorData = entoData.getNonVector();
    if (!nonVectorData.present())       // should be a validation error, but anyway...
      throw util::xml_scenario_error ("Neither vector nor non-vector data present in the XML!");
    if (util::ModelOptions::option( util::VECTOR_LIFE_CYCLE_MODEL ) ||
        util::ModelOptions::option( util::VECTOR_SIMPLE_MPD_MODEL ))
        throw util::xml_scenario_error("VECTOR_*_MODEL is only compatible with the vector model (and non-vector data is present).");
    model = new NonVectorModel(entoData, nonVectorData.get());
  }

  if( entoData.getScaledAnnualEIR().present() ){
      model->scaleEIR( entoData.getScaledAnnualEIR().get() / model->annualEIR );
      assert( vectors::approxEqual( model->annualEIR, entoData.getScaledAnnualEIR().get() ) );
  }

#ifdef WITHOUT_BOINC
  if( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ){
      //Note: after internal scaling (which doesn't imply exit)
      //but before external scaling.
      cout << "Total annual EIR: "<<model->annualEIR<<endl;
  }
#endif

  return model;
}


// The times here should be for the last updated index of arrays:
void TransmissionModel::ctsCbInputEIR (ostream& stream){
    //NOTE: because prevNow may be negative, we can't use mod_nn (hence moduloYearSteps):
    stream<<'\t'<<initialisationEIR[util::mod(sim::prevNow() / sim::oneTS(), sim::stepsPerYear())];
}
void TransmissionModel::ctsCbSimulatedEIR (ostream& stream){
    stream<<'\t'<<tsAdultEIR;
}
void TransmissionModel::ctsCbKappa (ostream& stream){
    // The latest time-step's kappa:
    stream<<'\t'<<laggedKappa[sim::now().moduloSteps(laggedKappa.size())];
}
void TransmissionModel::ctsCbNumTransmittingHumans (ostream& stream){
    stream<<'\t'<<numTransmittingHumans;
}


SimulationMode readMode(const string& str){
    if(str=="forced")
        return forcedEIR;
    else if(str=="dynamic")
        return dynamicEIR;
    else
        // Note: originally 3 (transientEIRknown) could be specified; now it's
        // set automatically.
        throw util::xml_scenario_error(string("mode attribute invalid: ").append(str));
}
TransmissionModel::TransmissionModel(const scnXml::Entomology& entoData) :
    simulationMode(forcedEIR),
    interventionMode(readMode(entoData.getMode())),
    laggedKappa(1, 0.0),        // if using non-vector model, it will resize this
    annualEIR(0.0),
    _annualAverageKappa(numeric_limits<double>::signaling_NaN()),
    _sumAnnualKappa(0.0),
    tsAdultEntoInocs(0.0),
    tsAdultEIR(0.0),
    surveyInputEIR(0.0),
    surveySimulatedEIR(0.0),
    adultAge(PerHost::adultAge()),
    numTransmittingHumans(0),
    tsNumAdults(0),
    timeStepNumEntoInocs (0)
{
  initialisationEIR.assign (sim::stepsPerYear(), 0.0);
  inoculationsPerAgeGroup.assign (Monitoring::AgeGroup::getNumGroups(), 0.0);
  timeStepEntoInocs.assign (Monitoring::AgeGroup::getNumGroups(), 0.0);

  // noOfAgeGroupsSharedMem must be at least as large as both of these to avoid
  // memory corruption or extra tests when setting/copying values
  noOfAgeGroupsSharedMem = std::max(Monitoring::AgeGroup::getNumGroups(), util::SharedGraphics::KappaArraySize);

  using Monitoring::Continuous;
  Continuous.registerCallback( "input EIR", "\tinput EIR", MakeDelegate( this, &TransmissionModel::ctsCbInputEIR ) );
  Continuous.registerCallback( "simulated EIR", "\tsimulated EIR", MakeDelegate( this, &TransmissionModel::ctsCbSimulatedEIR ) );
  Continuous.registerCallback( "human infectiousness", "\thuman infectiousness", MakeDelegate( this, &TransmissionModel::ctsCbKappa ) );
  Continuous.registerCallback( "num transmitting humans", "\tnum transmitting humans", MakeDelegate( this, &TransmissionModel::ctsCbNumTransmittingHumans ) );
}

TransmissionModel::~TransmissionModel () {
}


double TransmissionModel::updateKappa (const Population& population) {
    // We calculate kappa for output and non-vector model, and kappaByAge for
    // the shared graphics.

    double sumWt_kappa= 0.0;
    double sumWeight  = 0.0;
    kappaByAge.assign (noOfAgeGroupsSharedMem, 0.0);
    nByAge.assign (noOfAgeGroupsSharedMem, 0);
    numTransmittingHumans = 0;

    for (Population::ConstIter h = population.cbegin(); h != population.cend(); ++h) {
        //NOTE: calculate availability relative to age at end of time step;
        // not my preference but consistent with TransmissionModel::getEIR().
        double t = h->perHostTransmission.relativeAvailabilityHetAge(h->age(sim::ts1()).inYears());
        sumWeight += t;
        double pTransmit = 0.0;
        for( size_t g = 0; g < WithinHost::Genotypes::N(); ++g ){
            pTransmit += h->withinHostModel->probTransmissionToMosquito(
                h->getVaccine().getFactor( interventions::Vaccine::TBV ), g );
        }
        t *= pTransmit;
        sumWt_kappa += t;
        if( t > 0.0 )
            ++numTransmittingHumans;

        // kappaByAge and nByAge are used in the screensaver only
        Monitoring::AgeGroup ag = h->getMonitoringAgeGroup();
        kappaByAge[ag.i()] += t;
        ++nByAge[ag.i()];
    }


    size_t lKMod = sim::ts1().moduloSteps(laggedKappa.size());	// now
    if( population.size() == 0 ){     // this is valid
        laggedKappa[lKMod] = 0.0;        // no humans: no infectiousness
    } else {
        if ( !(sumWeight > DBL_MIN * 10.0) ){       // if approx. eq. 0, negative or an NaN
            ostringstream msg;
            msg<<"sumWeight is invalid: "<<sumWeight<<", "<<sumWt_kappa<<", "<<population.size();
            throw TRACED_EXCEPTION(msg.str(),util::Error::SumWeight);
        }
        laggedKappa[lKMod] = sumWt_kappa / sumWeight;
    }
    
    size_t tmod = sim::ts0().moduloYearSteps();
    
    //Calculate time-weighted average of kappa
    _sumAnnualKappa += laggedKappa[lKMod] * initialisationEIR[tmod];
    if (tmod == sim::stepsPerYear() - 1) {
        _annualAverageKappa = _sumAnnualKappa / annualEIR;	// inf or NaN when annualEIR is 0
        _sumAnnualKappa = 0.0;
    }

    // Shared graphics: report infectiousness
    if( mod_nn(sim::ts0(), sim::fromTS(6)) == sim::zero() ){
        for (size_t i = 0; i < noOfAgeGroupsSharedMem; i++)
            kappaByAge[i] /= nByAge[i];
        util::SharedGraphics::copyKappa(&kappaByAge[0]);
    }

    // Sum up inoculations this time step
    for (size_t group = 0; group < timeStepEntoInocs.size(); ++group) {
        inoculationsPerAgeGroup[group] += timeStepEntoInocs[group];
        // Reset to zero:
        timeStepEntoInocs[group] = 0.0;
    }
    timeStepNumEntoInocs = 0;

    tsAdultEIR = tsAdultEntoInocs / tsNumAdults;
    tsAdultEntoInocs = 0.0;
    tsNumAdults = 0;

    surveyInputEIR += initialisationEIR[tmod];
    surveySimulatedEIR += tsAdultEIR;
    
    return laggedKappa[lKMod];  // kappa now
}

double TransmissionModel::getEIR (Host::Human& human, SimTime age, double ageYears, OM::Monitoring::AgeGroup ageGroup) {
  /* For the NonVector model, the EIR should just be multiplied by the
   * availability. For the Vector model, the availability is also required
   * for internal calculations, but again the EIR should be multiplied by the
   * availability. */
  double EIR = calculateEIR (human, ageYears);

  //NOTE: timeStep*EntoInocs will rarely be used despite frequent updates here
  timeStepEntoInocs[ageGroup.i()] += EIR;
  timeStepNumEntoInocs ++;
  if( age >= adultAge ){
     tsAdultEntoInocs += EIR;
     tsNumAdults += 1;
  }
  util::streamValidate( EIR );
  return EIR;
}

void TransmissionModel::summarize () {
    Monitoring::Survey& survey = Monitoring::Survey::current();
    survey.setNumTransmittingHosts(laggedKappa[sim::now().moduloSteps(laggedKappa.size())]);
    survey.setAnnualAverageKappa(_annualAverageKappa);

    survey.setInoculationsPerAgeGroup (inoculationsPerAgeGroup);        // Array contents must be copied.
    inoculationsPerAgeGroup.assign (inoculationsPerAgeGroup.size(), 0.0);
    
    double duration = (sim::now() - lastSurveyTime).inSteps();
    if( duration == 0.0 ){
        if( !( surveyInputEIR == 0.0 && surveySimulatedEIR == 0.0 ) ){
            throw TRACED_EXCEPTION_DEFAULT( "non-zero EIR over zero duration??" );
        }
        duration = 1.0;   // avoid outputting NaNs. 0 isn't quite correct, but should do.
    }
    //TODO: we should really use bites-per-day or per year instead of per time
    // step. But we also can't just change an existing measure.
    survey.setInputEIR (surveyInputEIR / duration);
    survey.setSimulatedEIR (surveySimulatedEIR / duration);

    surveyInputEIR = 0.0;
    surveySimulatedEIR = 0.0;
    lastSurveyTime = sim::now();
}


// -----  checkpointing  -----

void TransmissionModel::checkpoint (istream& stream) {
    simulationMode & stream;
    interventionMode & stream;
    initialisationEIR & stream;
    laggedKappa & stream;
    annualEIR & stream;
    _annualAverageKappa & stream;
    _sumAnnualKappa & stream;
    adultAge & stream;
    tsAdultEntoInocs & stream;
    tsAdultEIR & stream;
    surveyInputEIR & stream;
    surveySimulatedEIR & stream;
    lastSurveyTime & stream;
    numTransmittingHumans & stream;
    tsNumAdults & stream;
    inoculationsPerAgeGroup & stream;
    timeStepEntoInocs & stream;
    timeStepNumEntoInocs & stream;
    noOfAgeGroupsSharedMem & stream;
    kappaByAge & stream;
    nByAge & stream;
}
void TransmissionModel::checkpoint (ostream& stream) {
    simulationMode & stream;
    interventionMode & stream;
    initialisationEIR & stream;
    laggedKappa & stream;
    annualEIR & stream;
    _annualAverageKappa & stream;
    _sumAnnualKappa & stream;
    adultAge & stream;
    tsAdultEntoInocs & stream;
    tsAdultEIR & stream;
    surveyInputEIR & stream;
    surveySimulatedEIR & stream;
    lastSurveyTime & stream;
    numTransmittingHumans & stream;
    tsNumAdults & stream;
    inoculationsPerAgeGroup & stream;
    timeStepEntoInocs & stream;
    timeStepNumEntoInocs & stream;
    noOfAgeGroupsSharedMem & stream;
    kappaByAge & stream;
    nByAge & stream;
}

} }

