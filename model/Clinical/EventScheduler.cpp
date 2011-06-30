/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#include "Clinical/EventScheduler.h"
#include "inputData.h"
#include "util/random.h"
#include "WithinHost/WithinHostModel.h"
#include "Transmission/PerHostTransmission.h"
#include "Monitoring/Surveys.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <limits>

namespace OM { namespace Clinical {
    using namespace OM::util;

TimeStep ClinicalEventScheduler::maxUCSeekingMemory(TimeStep::never);
TimeStep ClinicalEventScheduler::uncomplicatedCaseDuration(TimeStep::never);
TimeStep ClinicalEventScheduler::complicatedCaseDuration(TimeStep::never);
TimeStep ClinicalEventScheduler::extraDaysAtRisk(TimeStep::never);
vector<double> ClinicalEventScheduler::cumDailyPrImmUCTS;
double ClinicalEventScheduler::neg_v;

double ClinicalEventScheduler::hetWeightMultStdDev = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::minHetWeightMult = std::numeric_limits<double>::signaling_NaN();
AgeGroupInterpolation* ClinicalEventScheduler::weight = AgeGroupInterpolation::dummyObject();

double ClinicalEventScheduler::logOddsAbBase = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbNegTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbPosTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbNeed = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::logOddsAbInformal = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::oneMinusEfficacyAb = std::numeric_limits<double>::signaling_NaN();
AgeGroupInterpolation* ClinicalEventScheduler::severeNmfMortality = AgeGroupInterpolation::dummyObject();


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
    if (TimeStep::interval != 1)
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
    if (! (util::ModelOptions::option (util::INCLUDES_PK_PD)))
        throw util::xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");
    
    if( !InputData().getModel().getHuman().getWeight().present() ){
        throw util::xml_scenario_error( "model->human->weight element required by 1-day timestep model" );
    }
    weight = AgeGroupInterpolation::makeObject( InputData().getModel().getHuman().getWeight().get(), "weight" );
    hetWeightMultStdDev = InputData().getModel().getHuman().getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetWeightMult = 0.5 / weight->eval( 0.0 );
}

void ClinicalEventScheduler::setParameters (const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();
    
    maxUCSeekingMemory = TimeStep(coData.getMaxUCSeekingMemory());
    uncomplicatedCaseDuration = TimeStep(coData.getUncomplicatedCaseDuration());
    complicatedCaseDuration = TimeStep(coData.getComplicatedCaseDuration());
    extraDaysAtRisk = TimeStep(coData.getComplicatedRiskDuration()) - complicatedCaseDuration;
    if( uncomplicatedCaseDuration<TimeStep(1)
	|| complicatedCaseDuration<TimeStep(1)
	|| maxUCSeekingMemory<TimeStep(0)
	|| extraDaysAtRisk+complicatedCaseDuration<TimeStep(1)	// at risk at least 1 day
	|| extraDaysAtRisk>TimeStep(0)	// at risk longer than case duration
    ){
	throw util::xml_scenario_error(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }
    
    cumDailyPrImmUCTS.reserve( coData.getDailyPrImmUCTS().size() );
    double cumP = 0.0;
    for( scnXml::ClinicalOutcomes::DailyPrImmUCTSConstIterator it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it ){
        cumP += *it;
        cumDailyPrImmUCTS.push_back( cumP );
    }
    if( cumP < 0.999 || cumP > 1.001 ){
        throw util::xml_scenario_error( "Event scheduler: dailyPrImmUCTS seq must add up to 1" );
    }
    cumDailyPrImmUCTS.back() = 1.0;
    
    double alpha = exp( -InputData.getParameter( Params::CFR_NEG_LOG_ALPHA ) );
    if( !(0.0<=alpha && alpha<=1.0) ){
	throw util::xml_scenario_error(
            "Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    
    CaseManagementCommon::scaleCaseFatalityRate( alpha );
    neg_v = - InputData.getParameter( Params::CFR_SCALE_FACTOR );
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !esData.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();
        
        double pT = nmfDesc.getPrTreatment();
        logOddsAbBase = log( pT / (1.0 - pT) );
        logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        logOddsAbNeed = log(nmfDesc.getEffectNeed());
        logOddsAbInformal = log(nmfDesc.getEffectInformal());
        oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        severeNmfMortality = AgeGroupInterpolation::makeObject( nmfDesc.getCFR(), "CFR" );
    }
}

void ClinicalEventScheduler::cleanup () {
    AgeGroupInterpolation::freeObject( weight );
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
        ClinicalModel (cF),
        pgState (Pathogenesis::NONE),
        caseStartTime (TimeStep::never),
        timeOfRecovery (TimeStep::never),
        timeLastTreatment (TimeStep::never),
        previousDensity (numeric_limits<double>::quiet_NaN())
{
    if( tSF != 1.0 ){
	// p(treatment seeking) is part of tree & handled generically, so we
	// don't have a way of modifying it.
	throw xml_scenario_error("treatment seeking heterogeneity not supported");
    }
#ifndef NDEBUG
    int counter = 0;
#endif
    do {
        hetWeightMultiplier = util::random::gauss( 1.0, hetWeightMultStdDev );
#ifndef NDEBUG
        assert( counter < 100 );        // too many resamples: resamples should rarely be needed...
        ++counter;
#endif
    } while( hetWeightMultiplier < minHetWeightMult );
}
ClinicalEventScheduler::~ClinicalEventScheduler() {
    AgeGroupInterpolation::freeObject( severeNmfMortality );
}


// -----  other methods  -----

void ClinicalEventScheduler::massDrugAdministration(Human& human){
    // Note: we use the same medication method as with drugs as treatments, hence the actual
    // medication doesn't occur until the next timestep.
    // Note: we augment any existing medications, however future medications will replace any yet-
    // to-be-medicated MDA treatments (even all MDA doses when treatment happens immediately).
    ESCaseManagement::massDrugAdministration (
        ESHostData( human.getAgeInYears(), *human.withinHostModel, pgState ),
        medicateQueue, human.getInCohort(), human.getMonitoringAgeGroup()
    );
}

void ClinicalEventScheduler::doClinicalUpdate (Human& human, double ageYears){
    WithinHostModel& withinHostModel = *human.withinHostModel;
    // Run pathogenesisModel
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    util::streamValidate( (newState << 16) & pgState );
    
    if ( TimeStep::simulation1() == timeOfRecovery ) {
	if( pgState & Pathogenesis::DIRECT_DEATH ){
	    // Human dies this timestep (last day of risk of death)
	    _doomed = DOOMED_COMPLICATED;
	    
	    latestReport.update (human.getInCohort(), human.getMonitoringAgeGroup(), pgState);
        } else {
	    if ( pgState & Pathogenesis::COMPLICATED ) {
		if( random::uniform_01() < ESCaseManagement::pSequelaeInpatient( ageYears ) ){
		    pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
                }else{
		    pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
                }
	    } else {
		pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
            }
	    // report bout, at conclusion of episode:
	    latestReport.update (human.getInCohort(), human.getMonitoringAgeGroup(), pgState);
	    
	    // Individual recovers (and is immediately susceptible to new cases)
	    pgState = Pathogenesis::NONE;	// recovery (reset to healthy state)
	    
	    // And returns to transmission (if was removed)
	    human.perHostTransmission.removeFromTransmission( false );
	}
    }
    
    if ( newState & Pathogenesis::SICK ){
        // we have some new case: is it severe/complicated?
        if ( newState & Pathogenesis::COMPLICATED ){
            if ( pgState & Pathogenesis::COMPLICATED ) {
                // previously severe: no events happen for course of medication
            } else {
                // previously healthy or UC: progress to severe
                pgState = Pathogenesis::State (pgState | newState | Pathogenesis::RUN_CM_TREE);
                caseStartTime = TimeStep::simulation1();
            }
        } else {
            // uncomplicated case (UC/UC2/NMF): is it new?
            if (pgState & Pathogenesis::SICK) {
                // previously UC; nothing to do
            }else{
                if( caseStartTime < TimeStep::simulation1() ) {
                    // new UC case
                    pgState = Pathogenesis::State (pgState | newState | Pathogenesis::RUN_CM_TREE);
                    
                    double uVariate = random::uniform_01();
                    size_t i = 0;
                    for (; i < cumDailyPrImmUCTS.size(); ++i){
                        if( uVariate < cumDailyPrImmUCTS[i] ){
                            goto gotDelay;
                        }
                    }
                    assert(false);      // should have uVariate < 1 = cumDailyPrImmUCTS[len-1]
                    gotDelay:
                    // set start time: current time plus length of delay (days)
                    caseStartTime = TimeStep::simulation1() + TimeStep(i);
                }
            }
        }
        
        if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
            _doomed = -TimeStep::interval; // start indirect mortality countdown
    }
    
    if ( caseStartTime == TimeStep::simulation1() && pgState & Pathogenesis::RUN_CM_TREE ){
        // OK, we're about to run the CM tree
        pgState = Pathogenesis::State (pgState & ~Pathogenesis::RUN_CM_TREE);
        
	// If last treatment prescribed was in recent memory, consider second line.
	if (timeLastTreatment + Episode::healthSystemMemory > TimeStep::simulation1())
	    pgState = Pathogenesis::State (pgState | Pathogenesis::SECOND_CASE);
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
	CMAuxOutput auxOut = ESCaseManagement::execute(
	    ESHostData( ageYears, withinHostModel, pgState ), medicateQueue, human.getInCohort()
	);
	
        if( medicateQueue.size() ){	// I.E. some treatment was given
	    timeLastTreatment = TimeStep::simulation1();
            if( pgState & Pathogenesis::COMPLICATED ){
                Monitoring::Surveys.getSurvey(human.getInCohort())
                    .reportTreatments3( human.getMonitoringAgeGroup(), 1 );
            }else{
                if( pgState & Pathogenesis::SECOND_CASE ){
                    Monitoring::Surveys.getSurvey(human.getInCohort())
                        .reportTreatments2( human.getMonitoringAgeGroup(), 1 );
                }else{
                    Monitoring::Surveys.getSurvey(human.getInCohort())
                        .reportTreatments1( human.getMonitoringAgeGroup(), 1 );
                }
            }
        }
	
	if ( auxOut.hospitalisation != CMAuxOutput::NONE ) {	// in hospital
	    pgState = Pathogenesis::State (pgState | Pathogenesis::EVENT_IN_HOSPITAL);
	    
	    if( auxOut.hospitalisation == CMAuxOutput::DELAYED )
		++caseStartTime;
	}
	
	// Case fatality rate (first day of illness)
	// P(death) is some fixed input scaled by age-specific CFR.
	if( (pgState & Pathogenesis::COMPLICATED)
	    && !(pgState & Pathogenesis::DIRECT_DEATH)
	) {
	    double pDeath = CaseManagementCommon::caseFatality( ageYears );
	    // community fatality rate when not in hospital or delayed hospital entry
	    if( auxOut.hospitalisation != CMAuxOutput::IMMEDIATE )
		pDeath = CaseManagementCommon::getCommunityCaseFatalityRate( pDeath );
	    if (random::uniform_01() < pDeath) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH | Pathogenesis::EVENT_FIRST_DAY);
		// Human is killed at end of time at risk
		//timeOfRecovery += extraDaysAtRisk;	(no point updating; will be set later: ATORWD)
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
	
	if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
            if( (pgState & Pathogenesis::SICK) && !(pgState & Pathogenesis::COMPLICATED) ){
                // Have a NMF or UC malaria case
                double pNeedTreat = pathogenesisModel->pNmfRequiresTreatment( ageYears, (pgState & Pathogenesis::MALARIA) != Pathogenesis::NONE );
                bool needTreat = random::bernoulli(pNeedTreat);
                
                // Calculate chance of antibiotic administration:
                double pTreatment;
                if( auxOut.AB_provider == CMAuxOutput::NO_AB ){
                    // No treatment seeking
                    pTreatment = 0.0;
                }else{
                    // Treatment seeking; may be at a health facility or from
                    // the informal sector.
                    double logOddsAb = logOddsAbBase - logOddsAbNeed * pNeedTreat;
                    if( auxOut.AB_provider == CMAuxOutput::FACILITY ){
                        if( auxOut.diagnostic == CMAuxOutput::NEGATIVE ){
                            logOddsAb += logOddsAbNegTest;
                        }else if( auxOut.diagnostic == CMAuxOutput::POSITIVE ){
                            logOddsAb += logOddsAbPosTest;
                        }
                        if( needTreat ){
                            logOddsAb += logOddsAbNeed;
                        }
                    }else{
                        assert( auxOut.AB_provider == CMAuxOutput::INFORMAL );
                        logOddsAb += logOddsAbInformal;
                    }
                    
                    double oddsTreatment = exp( logOddsAb );
                    pTreatment = oddsTreatment / (1.0 + oddsTreatment);
                }
                
                double treatmentEffectMult = 1.0;
                
                if( random::uniform_01() < pTreatment ){
                    Monitoring::Surveys.getSurvey(human.getInCohort())
                        .reportAntibioticTreatments( human.getMonitoringAgeGroup(), 1 );
                    treatmentEffectMult = oneMinusEfficacyAb;
                }
                
                // In a severe NMF case (only when not malarial), there is a
                // chance of death:
                if( needTreat ){
                    double pDeath = severeNmfMortality->eval( ageYears ) * treatmentEffectMult;
                    if( random::uniform_01() < pDeath ){
                        pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
                    }
                }
            }
        }
    } else {
	// No new event (haven't changed state this timestep).
	
	// Case fatality rate (subsequent days)
	// Complicated case & at risk of death (note: extraDaysAtRisk <= 0)
	if( (pgState & Pathogenesis::COMPLICATED)
	    && !(pgState & Pathogenesis::DIRECT_DEATH)
	    && (TimeStep::simulation1() < timeOfRecovery + extraDaysAtRisk)
	) {
	    // In complicated episodes, S(t), the probability of survival on
	    // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
	    // for parasite density Y(t). v_neg below is -v.
	    if( withinHostModel.getTotalDensity() > 0.0 ){	// avoid division by zero
		double parasiteReductionEffect = withinHostModel.getTotalDensity() / previousDensity;
		double pDeath = 1.0 - exp( neg_v * parasiteReductionEffect );
		// community fatality rate when not in hospital
		if( !(pgState & Pathogenesis::EVENT_IN_HOSPITAL) )
		    pDeath = CaseManagementCommon::getCommunityCaseFatalityRate( pDeath );
		if (random::uniform_01() < pDeath) {
		    pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		    // Human is killed at end of time at risk
		    timeOfRecovery += extraDaysAtRisk;	// may be re-set later (see ATORWD)
		}
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
    }
    
    // Start of case. Not necessarily start of sickness due to treatment-seeking
    // delays and travel time.
    if( caseStartTime == TimeStep::simulation1() ){
	// Patients in hospital are removed from the transmission cycle.
	// This should have an effect from the start of the next timestep.
	// NOTE: This is not very accurate, but considered of little importance.
	if (pgState & Pathogenesis::EVENT_IN_HOSPITAL)
	    human.perHostTransmission.removeFromTransmission( true );
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    // complicatedCaseDuration should to some respects be associated
	    // with medication duration, however ongoing medications after
	    // exiting hospital are OK and medications terminating before the
	    // end of hospitalisation shouldn't matter too much if the person
	    // can't recieve new infections due to zero transmission in hospital.
	    timeOfRecovery = TimeStep::simulation1() + complicatedCaseDuration;
	    // Time should be adjusted to end of at-risk period when patient dies:
	    if( pgState & Pathogenesis::DIRECT_DEATH )	// death may already have been determined
		timeOfRecovery += extraDaysAtRisk;	// ATORWD (search keyword)
	} else {
	    timeOfRecovery = TimeStep::simulation1() + uncomplicatedCaseDuration;
	}
    }
    
    
    // Process pending medications (in interal queue) and apply/update:
    for (list<MedicateData>::iterator it = medicateQueue.begin(); it != medicateQueue.end();) {
	list<MedicateData>::iterator next = it;
	++next;
        if ( it->time < 1.0 ) { // Medicate today's medications
            double bodyMass = ageToWeight( ageYears );
	    withinHostModel.medicate (it->abbrev, it->qty, it->time, it->duration, bodyMass);
            if( it->duration > 0.0 ){
                Monitoring::Surveys.getSurvey(human.getInCohort())
                .report_Clinical_DrugUsageIV (it->abbrev, it->cost_qty * bodyMass);
            }else{      // 0 or NaN
                Monitoring::Surveys.getSurvey(human.getInCohort())
                .report_Clinical_DrugUsage (it->abbrev, it->cost_qty);
            }
	    medicateQueue.erase (it);
        } else {   // and decrement treatment seeking delay for the rest
            it->time -= 1.0;
        }
	it = next;
    }
    
    if( human.cohortFirstTreatmentOnly && timeLastTreatment == TimeStep::simulation1() ){
        human.removeFromCohort();
    }
    if( human.cohortFirstBoutOnly && (pgState & Pathogenesis::SICK) ){
        human.removeFromCohort();
    }
}


void ClinicalEventScheduler::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    int s;
    s & stream;
    pgState = Pathogenesis::State(s);
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
    hetWeightMultiplier & stream;
    medicateQueue & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
    hetWeightMultiplier & stream;
    medicateQueue & stream;
}

} }