/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "mon/info.h"
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

  if( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ){
      //Note: after internal scaling (which doesn't imply exit)
      //but before external scaling.
      cout << "Total annual EIR: "<<model->annualEIR<<endl;
  }

  return model;
}


// The times here should be for the last updated index of arrays:
void TransmissionModel::ctsCbInputEIR (ostream& stream){
    int prevStep = (sim::now() - SimTime::oneTS()) / SimTime::oneTS();
    //Note: prevStep may be negative, hence util::mod not mod_nn:
    stream<<'\t'<<initialisationEIR[util::mod(prevStep, SimTime::stepsPerYear())];
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
TransmissionModel::TransmissionModel(const scnXml::Entomology& entoData,
                                     size_t nGenotypes) :
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
    tsNumAdults(0)
{
    initialisationEIR.assign (SimTime::stepsPerYear(), 0.0);
    
  using Monitoring::Continuous;
  Continuous.registerCallback( "input EIR", "\tinput EIR", MakeDelegate( this, &TransmissionModel::ctsCbInputEIR ) );
  Continuous.registerCallback( "simulated EIR", "\tsimulated EIR", MakeDelegate( this, &TransmissionModel::ctsCbSimulatedEIR ) );
  Continuous.registerCallback( "human infectiousness", "\thuman infectiousness", MakeDelegate( this, &TransmissionModel::ctsCbKappa ) );
  Continuous.registerCallback( "num transmitting humans", "\tnum transmitting humans", MakeDelegate( this, &TransmissionModel::ctsCbNumTransmittingHumans ) );
}

TransmissionModel::~TransmissionModel () {
}


double TransmissionModel::updateKappa () {
    // We calculate kappa for output and the non-vector model.
    double sumWt_kappa= 0.0;
    double sumWeight  = 0.0;
    numTransmittingHumans = 0;

    foreach(const Host::Human& human, sim::humanPop().crange()) {
        //NOTE: calculate availability relative to age at end of time step;
        // not my preference but consistent with TransmissionModel::getEIR().
        const double avail = human.perHostTransmission.relativeAvailabilityHetAge(
            human.age(sim::ts1()).inYears());
        sumWeight += avail;
        const double tbvFactor = human.getVaccine().getFactor( interventions::Vaccine::TBV );
        const double pTransmit = human.withinHostModel->probTransmissionToMosquito( tbvFactor, 0 );
        const double riskTrans = avail * pTransmit;
        sumWt_kappa += riskTrans;
        if( riskTrans > 0.0 )
            ++numTransmittingHumans;
    }


    size_t lKMod = sim::ts1().moduloSteps(laggedKappa.size());	// now
    if( sim::humanPop().size() == 0 ){     // this is valid
        laggedKappa[lKMod] = 0.0;        // no humans: no infectiousness
    } else {
        if ( !(sumWeight > DBL_MIN * 10.0) ){       // if approx. eq. 0, negative or an NaN
            ostringstream msg;
            msg<<"sumWeight is invalid: "<<sumWeight<<", "<<sumWt_kappa
                    <<", "<<sim::humanPop().size();
            throw TRACED_EXCEPTION(msg.str(),util::Error::SumWeight);
        }
        laggedKappa[lKMod] = sumWt_kappa / sumWeight;
    }
    
    size_t tmod = sim::ts0().moduloYearSteps();
    
    //Calculate time-weighted average of kappa
    _sumAnnualKappa += laggedKappa[lKMod] * initialisationEIR[tmod];
    if (tmod == SimTime::stepsPerYear() - 1) {
        _annualAverageKappa = _sumAnnualKappa / annualEIR;	// inf or NaN when annualEIR is 0
        _sumAnnualKappa = 0.0;
    }
    
    tsAdultEIR = tsAdultEntoInocs / tsNumAdults;
    tsAdultEntoInocs = 0.0;
    tsNumAdults = 0;
    
    surveyInputEIR += initialisationEIR[tmod];
    surveySimulatedEIR += tsAdultEIR;
    
    return laggedKappa[lKMod];  // kappa now
}

double TransmissionModel::getEIR( Host::Human& human, SimTime age,
                    double ageYears, vector<double>& EIR )
{
    /* For the NonVector model, the EIR should just be multiplied by the
     * availability. For the Vector model, the availability is also required
     * for internal calculations, but again the EIR should be multiplied by the
     * availability. */
    calculateEIR( human, ageYears, EIR );
    util::streamValidate( EIR );
    
    double allEIR = vectors::sum( EIR );
    if( age >= adultAge ){
        tsAdultEntoInocs += allEIR;
        tsNumAdults += 1;
    }
    return allEIR;
}

void TransmissionModel::summarize () {
    mon::reportStatMF( mon::MVF_NUM_TRANSMIT, laggedKappa[sim::now().moduloSteps(laggedKappa.size())] );
    mon::reportStatMF( mon::MVF_ANN_AVG_K, _annualAverageKappa );
    
    if( !mon::isReported() ) return;    // cannot use counters below when not reporting
    
    double duration = (sim::now() - lastSurveyTime).inSteps();
    if( duration == 0.0 ){
        assert( surveyInputEIR == 0.0 && surveySimulatedEIR == 0.0 );
        duration = 1.0;   // avoid outputting NaNs. 0 isn't quite correct, but should do.
    }
    mon::reportStatMF( mon::MVF_INPUT_EIR, surveyInputEIR / duration );
    mon::reportStatMF( mon::MVF_SIM_EIR, surveySimulatedEIR / duration );
    
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
}

} }

