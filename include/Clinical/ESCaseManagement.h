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

#ifndef Hmod_ESCaseManagement
#define Hmod_ESCaseManagement

#include "Global.h"
#include "Clinical/CMDecisionTree.h"    // needed for ESDecisionMap
#include "WithinHost/WHInterface.h"
#include "Monitoring/Survey.h"
#include "schema/pharmacology.h"

#include <cassert>
#include <list>
#include <limits>
#include <boost/ptr_container/ptr_unordered_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class ESCaseManagementSuite;
class ESDecisionTreeSuite;

namespace scnXml{
    
}
namespace OM { namespace Clinical {
using WithinHost::WHInterface;
using Monitoring::Survey;


/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *****************************************************************************/
class ESCaseManagement {
public:
    /** Load health system data from initial data or an intervention's data (both from XML).
    * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem (const scnXml::HealthSystem& healthSystem);
    
    /** Set up MDA drug. Must be called if massDrugAdministration() is
        * ever used to deploy an MDA intervention. */
    static void initMDA (const scnXml::DecisionTree& desc);
    
    static void massDrugAdministration(
        const CMHostData& hostData,
        const Host::Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport
    );
    
    static CMDTOut execute( const CMHostData& hostData );
    
private:
    static auto_ptr<CMDecisionTree> uncomplicated, complicated, mda;
};

} }
#endif