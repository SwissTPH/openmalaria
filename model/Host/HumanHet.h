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

// This is a simple header, for direct inclusion in Human.cpp only.

#include "util/random.h"
#include "util/ModelOptions.h"

namespace OM {
namespace Host {

/// (Some) heterogeneity parameters of humans
struct HumanHet {
    /* Human heterogeneity; affects:
     * comorbidityFactor (stored in PathogenesisModel)
     * treatmentSeekingFactor (stored in CaseManagementModel)
     * availabilityFactor (stored in Transmission::PerHost) */
    double comorbidityFactor;
    double treatmentSeekingFactor;
    double availabilityFactor;
    
    HumanHet() : comorbidityFactor(1.0), treatmentSeekingFactor(1.0),
        availabilityFactor(1.0) {}
    
    static HumanHet sample(util::LocalRng& rng)
    {
        HumanHet het;
        if(util::ModelOptions::option (util::TRANS_HET)){
            het.availabilityFactor = 0.2;
            if( rng.bernoulli(0.5) ){
                het.availabilityFactor = 1.8;
            }
        }
        if(util::ModelOptions::option (util::COMORB_HET)){
            het.comorbidityFactor = 0.2;
            if( rng.bernoulli(0.5) ){
                het.comorbidityFactor = 1.8;
            }
        }
        if(util::ModelOptions::option (util::TREAT_HET)){
            het.treatmentSeekingFactor = 0.2;
            if( rng.bernoulli(0.5) ){
                het.treatmentSeekingFactor = 1.8;
            }
        }
        if(util::ModelOptions::option (util::TRANS_TREAT_HET)){
            het.treatmentSeekingFactor = 0.2;
            het.availabilityFactor = 1.8;
            if( rng.bernoulli(0.5) ){
                het.treatmentSeekingFactor = 1.8;
                het.availabilityFactor = 0.2;
            }
        }else if(util::ModelOptions::option (util::COMORB_TREAT_HET)){
            if( rng.bernoulli(0.5) ){
                het.comorbidityFactor = 1.8;
                het.treatmentSeekingFactor = 0.2;
            }else{
                het.comorbidityFactor = 0.2;
                het.treatmentSeekingFactor = 1.8;
            }
        }else if(util::ModelOptions::option (util::COMORB_TRANS_HET)){
            het.availabilityFactor = 1.8;
            het.comorbidityFactor = 1.8;
            if( rng.bernoulli(0.5) ){
                het.availabilityFactor = 0.2;
                het.comorbidityFactor = 0.2;
            }
        }else if(util::ModelOptions::option (util::TRIPLE_HET)){
            het.availabilityFactor = 1.8;
            het.comorbidityFactor = 1.8;
            het.treatmentSeekingFactor = 0.2;
            if( rng.bernoulli(0.5) ){
                het.availabilityFactor = 0.2;
                het.comorbidityFactor = 0.2;
                het.treatmentSeekingFactor = 1.8;
            }
        }
        return het;
    }

};

}
}
