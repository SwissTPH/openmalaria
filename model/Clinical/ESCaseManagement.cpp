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

#include "Clinical/ESCaseManagement.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/EventScheduler.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"
#include "util/ModelOptions.h"

#include <set>
#include <sstream>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

namespace OM { namespace Clinical {
    using namespace OM::util;
    using namespace boost::assign;
    using namespace Monitoring;
    using boost::format;

// -----  ESTreatmentSchedule  -----

ESTreatmentSchedule::ESTreatmentSchedule (const scnXml::PKPDSchedule& sched) {
    const ::scnXml::PKPDSchedule::MedicateSequence& mSeq = sched.getMedicate();
    medications.resize (mSeq.size ());
    for (size_t j = 0; j < mSeq.size(); ++j) {
	medications[j].abbrev = mSeq[j].getDrug();
	medications[j].qty = mSeq[j].getMg();
	medications[j].cost_qty = medications[j].qty;	// effective quantity w.r.t. cost starts off equal to quantity
	medications[j].time = mSeq[j].getHour() / 24.0;	// convert from hours to days
	if( mSeq[j].getDuration().present() ){
	    if( !(mSeq[j].getDuration().get() > 0.0) ){
		throw util::xml_scenario_error( "duration of an IV dose must be some positive amount of time" );
	    }
	    medications[j].duration = mSeq[j].getDuration().get() / 24.0;
	}
    }
}


// -----  ESDecisionMap  -----

void ESDecisionMap::initialize (const ::scnXml::DecisionTree& xmlCM, TreeType treeType, bool reinitialise) {
    //TODO: set up decisions
}

ESDecisionMap::~ESDecisionMap () {
}

void ESDecisionMap::execute( CMHostData hostData ) const{
    //TODO
}


// -----  ESCaseManagement  -----

ESDecisionMap ESCaseManagement::uncomplicated,
    ESCaseManagement::complicated;
ESDecisionMap ESCaseManagement::mda;
// ESDecisionMap decisionsUC, decisionsSev, decisionsMDA;

void ESCaseManagement::setHealthSystem(
    const scnXml::HealthSystem& healthSystem)
{
    if( !healthSystem.getEventScheduler().present() )
	throw util::xml_scenario_error ("Expected EventScheduler section in "
        "healthSystem data (initial or intervention)");
    const scnXml::HSEventScheduler& esData =
        healthSystem.getEventScheduler().get();
    uncomplicated.initialize (esData.getUncomplicated (),
                              ESDecisionMap::Uncomplicated, true);
    complicated.initialize (esData.getComplicated (),
                            ESDecisionMap::Complicated, true);
    
    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters( esData );
}

void ESCaseManagement::initMDA (const scnXml::DecisionTree& desc){
    mda.initialize( desc, ESDecisionMap::MDA, false );
}

void ESCaseManagement::massDrugAdministration(
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue,
        const Host::Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport
){
    Survey::current().addInt( screeningReport, human, 1 );
    mda.execute( hostData );
    if( true /*FIXME anyTreatment*/ ){
        Survey::current().addInt( drugReport, human, 1 );
    }
}

CMAuxOutput ESCaseManagement::execute (
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue
) {
    assert (hostData.pgState & Episode::SICK);
    // We always remove any queued medications.
    medicateQueue.clear();
    
    ESDecisionMap *map = (hostData.pgState & Episode::COMPLICATED) ?
        &complicated : &uncomplicated;
    map->execute( hostData );
    
    CMAuxOutput auxOut;
    auxOut.hospitalisation = CMAuxOutput::NONE;
    /*FIXME
    if( hostData.pgState & Episode::COMPLICATED )
	auxOut.hospitalisation = map->hospitalisation(outcome);
    auxOut.diagnostic = map->diagnostic(outcome);
    auxOut.AB_provider = map->AB_provider(outcome);
    */
    return auxOut;
}

} }