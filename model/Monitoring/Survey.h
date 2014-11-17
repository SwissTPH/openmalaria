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

#ifndef Hmod_Survey
#define Hmod_Survey

#include "Global.h"
#include "Monitoring/SurveyMeasure.h"
#include "Monitoring/AgeGroup.h"        // only needed for special version of addInt() used by Episode
#include "WithinHost/Diagnostic.h"
#include "Parameters.h"
#include "util/checkpoint_containers.h"
#include <bitset>
#include <map>

class UnittestUtil;
namespace scnXml{ class Monitoring; }
namespace OM {
namespace Host {
    class Human;
}
namespace Monitoring {
using WithinHost::Diagnostic;

/// Data struct for a single survey.
class Survey {
private:
    
    ///@brief Static members (options from XML). Parameters only set by init().
    //@{
    /// Initialize static parameters.
    static void init( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring,
                   size_t nSurveys );
    
    /// Encoding of which summary options are active in XML is converted into
    /// this array for easier reading (and to make changing encoding within XML easier).
    static bitset<SM::NUM_SURVEY_OPTIONS> active;
    //@}
  
public:
    // Constructor used by SurveysType. Call allocate() explicitly for allocation.
    Survey();
    
    ///@brief Static access functions: these are here so that most users don't need to include Surveys.h
    //@{
    /** Get access to the current survey.
     * 
     * This is an inline function and should be very fast to look up.
     * 
     * Note: current() == getSurvey(getSurveyNumber()) */
    inline static Survey& current(){ return *m_current; }
    
    /** Returns the number of the current survey. Use this to report
     * retrospectively. */
    inline static size_t getSurveyNumber(){ return m_surveyNumber; }
    
    /** Return Survey number n (counting from 1). Use this along with
     * getSurveyNumber() to report retrospectively; in most cases this is not
     * needed and current() can be used instead. */
    static Survey& getSurvey(size_t n);
    
    /** Return the time of the final survey.
     *
     * We use this to control when the simulation ends.
     * This isn't quite the same as before when the simulation end was
     * explicitly specified and has a small affect on
     * infantAllCauseMortality (survey 21) output. */
    static SimTime getLastSurveyTime ();
    
    static inline const Diagnostic& diagnostic(){
        assert( m_diagnostic != 0 );
        return *m_diagnostic;
    }
    
    /** Humans should store a "cohort set" identifier which is initially 0.
     * Whenever a human gains or loses membership status in some
     * sup-population, it should update that value with this function.
     * 
     * @param old       Old identifier value (initially 0)
     * @param subPop    Sub-population to which membership status changed
     * @param isMember  New membership status
     * @returns         New identifier value
     */
    static uint32_t updateCohortSet( uint32_t old,
        interventions::ComponentId subPop, bool isMember );
    //@}
    
    ///@brief Set outputs without extra categorisation
    //@{
    /** Number of hosts transmitting to mosquitoes, reported as nTransmit. */
    void setNumTransmittingHosts(double value) { m_nTransmit = value; }
    /** Reported as annAvgK **/
    void setAnnualAverageKappa(double kappa) { m_annAvgK = kappa; }
    void setInputEIR (double v) { m_inputEIR = v; }
    void setSimulatedEIR (double v) { m_simulatedEIR = v; }
    //@}
    
    ///@brief Set outputs per vector species
    //@{
    Survey& set_Vector_Nv0 (string key, double v) { data_Vector_Nv0[key] = v; return *this; }
    Survey& set_Vector_Nv (string key, double v) { data_Vector_Nv[key] = v; return *this; }
    Survey& set_Vector_Ov (string key, double v) { data_Vector_Ov[key] = v; return *this; }
    Survey& set_Vector_Sv (string key, double v) { data_Vector_Sv[key] = v; return *this; }
    //@}
    
    void setInoculationsPerAgeGroup (vector<double>& v) {
        m_inoculationsPerAgeGroup = v;	// copies v, not just its reference
    }
  
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        m_nTransmit & stream;
        m_annAvgK & stream;
        m_inputEIR & stream;
        m_simulatedEIR & stream;
        
        data_Vector_Nv0 & stream;
        data_Vector_Nv & stream;
        data_Vector_Ov & stream;
        data_Vector_Sv & stream;
        
        m_inoculationsPerAgeGroup & stream;
  }
  
private:
    /** Resizes all vectors, allocating memory.
     * 
     * This is a separate initialisation step to make allocation explicit and
     * avoid accidental allocations when manipulating containers of Survey
     * elements. */
    void allocate ();
    
    /** Write out arrays
     * @param outputFile Stream to write to
     * @param survey Survey number (starting from 1) */
    void writeSummaryArrays (ostream& outputFile, int survey);
    
    /// @brief Data stored for reporting; all of this is per survey
    //@{
    // no further categorisation:
    double m_nTransmit;
    double m_annAvgK;
    double m_inputEIR;
    double m_simulatedEIR;
    
    // data categorised by vector species:
    map<string,double> data_Vector_Nv0;
    map<string,double> data_Vector_Nv;
    map<string,double> data_Vector_Ov;
    map<string,double> data_Vector_Sv;
    
    // data categorised by human age group:
    vector<double> m_inoculationsPerAgeGroup;
    //@}
    
    // ———  static members  ———
    
    /** Index for the time dimention of the summary arrays
     * Index starts from 1 for used surveys; is 0 to write to dummy survey. */
    static size_t m_surveyNumber;
    
    /** Points to SurveysType::surveys[m_surveyNumber] (or the dummy element 
     * SurveysType::surveys[0] before the intervention period and after
     * completion of last survey). This is for data being collected for the
     * next survey. */
    static Survey *m_current;
    
    static const Diagnostic* m_diagnostic;
    
    friend class SurveysType;
    friend class ::UnittestUtil;
};

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

} }
#endif
