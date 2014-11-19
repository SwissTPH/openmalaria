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

#ifndef Hmod_Surveys
#define Hmod_Surveys

#include "Monitoring/Survey.h"
#include "Global.h"
#include "Host/Human.h"
#include <assert.h>

namespace OM {
namespace interventions {
    class InterventionManager;
}
namespace Monitoring {
    
/** Class to collect surveys and write them out.
 *
 * Surveys are written to the file output.txt. There is a one-to-one mapping of
 * lines to data entries (except the file ends with a new line). Data columns are
 * separated by tabs on each.
 * 
 * The first column lists the survey number as an integer, counting from one,
 * the second column a "group" parameter as a string (precise meaning depends
 * on the measure), the third column the measure as an integer ID (the values
 * in the SurveyCodes enum), and the forth a value (integer or floating-point,
 * but when exported to the database always considered a double). */
class SurveysType
{
public:
    ///@brief Init, output, checkpointing functions
    //@{
    /** Second initialisation step: must happen after the InterventionManager is set up. */
    void init2( const scnXml::Monitoring& monitoring );
    
    //! Write all the summary arrays requested by summaryOption to output.txt
    void writeSummaryArrays();
    //@}
    
    ///@brief Simple getters
    //@{
    /** Get the number of cohort sets (i.e. two to the power of the number of
     * sub-populations considered cohorts). */
    uint32_t numCohortSets()const;
    
    /** Get the output cohort set numeric identifier given the internal one
     * (as returned by Survey::updateCohortSet()). */
    uint32_t cohortSetOutputId( uint32_t cohortSet )const;
    //@}
    
private:
    friend class Survey;
};
/// Data — entry-point for using Surveys. Checkpointed.
extern SurveysType Surveys;
} }
#endif
