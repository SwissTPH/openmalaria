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

#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CaseManagementCommon.h"
#include "Clinical/CMDecisionTree.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/WHVivax.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;
using namespace Monitoring;

enum CaseType { FirstLine, SecondLine, NumCaseTypes };

// These parameters are set by setHealthSystem() and do not need checkpointing.
auto_ptr<CMDecisionTree> treeUCOfficial, treeUCSelfTreat;

// ———  static, init  ———

void DecisionTree5Day::setHealthSystem(const scnXml::HSDT5Day& hsDescription){
    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    accessUCSelfTreat[SecondLine] = 0.0;
    accessUCAny[FirstLine] = accessUCOfficial1 + accessUCSelfTreat[FirstLine];
    accessUCAny[SecondLine] = accessUCOfficial1 + accessUCSelfTreat[SecondLine];
    accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();
    
    if( accessUCOfficial1 < 0.0 || accessUCSelfTreat[FirstLine] < 0.0 || accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || accessUCSelfTreat[SecondLine] < 0.0 || accessUCAny[SecondLine] > 1.0 ||
        accessSevere < 0.0 || accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }
    
    treeUCOfficial.reset( CMDecisionTree::create( hsDescription.getTreeUCOfficial(), true ).release() );
    treeUCSelfTreat.reset( CMDecisionTree::create( hsDescription.getTreeUCSelfTreat(), true ).release() );
    
    cureRateSevere = hsDescription.getCureRateSevere().getValue();
    treatmentSevere = WHInterface::addTreatment( hsDescription.getTreatmentSevere() );
    
    if( hsDescription.getPrimaquine().present() ){
        if( !ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
            throw util::xml_scenario_error( "health-system's primaquine element only supported by vivax" );
        WithinHost::WHVivax::setHSParameters( hsDescription.getPrimaquine().get() );
    }
}


// ———  per-human, update  ———

void DecisionTree5Day::uncomplicatedEvent ( Human& human, Episode::State pgState ){
    latestReport.update (human, Episode::State( pgState ) );
    
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (m_tLastTreatment + healthSystemMemory > sim::now1()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }
    
    double x = random::uniform_01();
    if( x < accessUCAny[regimen] * m_treatmentSeekingFactor ){
        CMHostData hostData( human, human.getAgeInYears(), pgState );
        
        // Run tree (which may deploy treatment)
        CMDTOut output = ( x < accessUCSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            treeUCSelfTreat->exec( hostData ) :
            treeUCOfficial->exec( hostData );
        
        if( output.treated ){   // if any treatment or intervention deployed
            m_tLastTreatment = sim::now1();
            Survey::current().addInt( measures[regimen], human, 1 );
        }
        
        if( human.withinHostModel->optionalPqTreatment() )
            Survey::current().addInt( Report::MI_PQ_TREATMENTS, human, 1 );
    } else {
        // No care sought
    }
}

}
}
