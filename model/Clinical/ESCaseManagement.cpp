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


// -----  ESCaseManagement  -----

auto_ptr<CMDecisionTree> ESCaseManagement::uncomplicated,
    ESCaseManagement::complicated,
    ESCaseManagement::mda;

void ESCaseManagement::setHealthSystem(const scnXml::HSEventScheduler& esData){
    uncomplicated = CMDecisionTree::create( esData.getUncomplicated () );
    complicated = CMDecisionTree::create( esData.getComplicated() );
    
    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters( esData );
}

void ESCaseManagement::initMDA (const scnXml::DecisionTree& desc){
    mda = CMDecisionTree::create( desc );
}

void ESCaseManagement::massDrugAdministration(
        const CMHostData& hostData,
        const Host::Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport
){
    Survey::current().addInt( screeningReport, human, 1 );
    CMDTOut out = mda->exec( hostData );
    if( out.treated ){
        Survey::current().addInt( drugReport, human, 1 );
    }
}

CMDTOut ESCaseManagement::execute (
        const CMHostData& hostData
) {
    assert (hostData.pgState & Episode::SICK);
    //TODO: should we remove any existing prescriptions?
    //TODO: Note that these trees do both "access" and "case management" decisions.
    //TODO: medicateQueue.clear();
    
    CMDecisionTree& tree = (hostData.pgState & Episode::COMPLICATED) ?
        *complicated : *uncomplicated;
    return tree.exec( hostData );
}

} }