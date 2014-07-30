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

#ifndef Hmod_PathogenesisModel
#define Hmod_PathogenesisModel

#include "Global.h"
#include "Parameters.h"
#include "Monitoring/Survey.h"
#include "WithinHost/Pathogenesis/State.h"
#include "util/AgeGroupInterpolation.h"

namespace scnXml{
    class HSESNMF;
    class Clinical;
}

namespace OM { namespace WithinHost { namespace Pathogenesis {
using util::AgeGroupInterpolator;

/*! PathogenesisModel abstract base class.
 *
 * Previously named MorbidityModel and PresentationModel. */
class PathogenesisModel {
public:
    /// Calls static init on correct PathogenesisModel.
    static void init( const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly );

    /** Create a sub-class instance, dependant on global options.
    *
    * @param cF = Comorbidity factor (currently set in Human). */
    static PathogenesisModel* createPathogenesisModel(double cF);

    
    virtual ~PathogenesisModel() {}

    /** Determines the health of the individual based on his/her parasitemia.
     *
     * May introduce severe or uncomplicated cases of malaria, as well as non-
     * malaria fevers. */
    StatePair determineState(double ageYears, double timeStepMaxDensity, double endDensity);
    
    /** For Vivax: determine the chance of a NMF and sample, returning either
     * NONE or STATE_NMF. */
    static Pathogenesis::State sampleNMF( double ageYears );
    
    /** Summarize PathogenesisModel details
     *
     * Only PyrogenPathogenesis implements this; other models don't have anything
     * to add to the summary. */
    virtual void summarize (const Host::Human& human) {}

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }

protected:
    /** Create a PathogenesisModel. */
    PathogenesisModel(double cF);

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

private:
    ///@brief static vars set by init()
    //@{
    /// Comorbidity prevalence at birth as a risk factor for indirect mortality
    static double indirRiskCoFactor_18;
    /// sevMal: critical density for severe malaria bout (Y*B1)
    static double sevMal_21;
    /// Critical age for co-morbidity (for both severe and indirect)
    static double critAgeComorb_30;
    /// Comorbidity prevalence at birth as a risk factor for severe
    static double comorbintercept_24;
    /// Rate of Non-Malaria Fever incidence by age. Non-seasonal.
    static AgeGroupInterpolator NMF_incidence;
    //@}

protected:
    virtual double getPEpisode(double timeStepMaxDensity, double totalDensity)=0;

    //! comorbidity factor for heterogeneity
    double _comorbidityFactor;
};

} } }
#endif
