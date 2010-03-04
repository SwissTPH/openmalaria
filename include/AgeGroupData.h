/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_AgeGroupData
#define Hmod_AgeGroupData

#include <cstddef>	// for size_t
#include <map>
#include <vector>

using namespace std;

/** Encapsulation for data according to some old age-groups. */
class AgeGroupData
{
public:
    /// Initialization
    static void initParameters ();
    
    static double getAgeSpecificRelativeAvailability (double ageYears);
    static double ageToWeight (double ageYears);
    
private:
    /// Get the appropriate index within agemax, ageSpecificRelativeAvailability
    /// and wtprop for this age (in years). Also used by PerHostTransmission.
    static size_t getAgeGroup (double age);
    static const map<double,size_t> fillAgeGroups();

    ///@brief Age-group variables for wtprop and ageSpecificRelativeAvailability
    //@{
    //! Number of age groups to use
    static const size_t nages= 22;
    //! Maximum of each age category
    static const double agemax[nages];

    static const map<double,size_t> ageMap;
  
    /** Average number of bites for each age as a proportion of the maximum.
    *
    * Set by constructor from constants (bsa_prop). */
    static double ageSpecificRelativeAvailability[nages];

    //! Proportionate body surface area
    /* 
    The body surface area is expressed as proportions of 0.5*those in 
    the reference age group.In some models we have used calculations of weight and in others surface area, based on 
    Mosteller RD: Simplified Calculation of Body Surface Area. N Engl J Med 1987 Oct 22;317(17):1098 (letter) 
    These values are retained here should they be required for future comparisons 
    */ 
    static const double bsa_prop[nages];
    
    //! Relative weights by age group
    /** Relative weights, based on data in InputTables\wt_bites.csv 
    The data are for Kilombero, Tanzania, taken from the Keiser et al (diploma
    thesis). The original source was anthropometric studies by Inez Azevedo Reads
    in weights by age group. The weights are expressed as proportions of 0.5*those
    in the reference age group. */
    static const double wtprop[nages];
    //@}
};

#endif
