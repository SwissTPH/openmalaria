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

#include "AgeGroupData.h"
#include "Global.h"

// NOTE: I'd rather use these than x.99 values, but it changes things!
//const double AgeGroupData::agemin[nages] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60 };
const double AgeGroupData::agemax[nages] = { 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 19.99, 24.99, 29.99, 39.99, 49.99, 59.99, 60.99 };

const double AgeGroupData::bsa_prop[AgeGroupData::nages] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
double AgeGroupData::ageSpecificRelativeAvailability[AgeGroupData::nages];

const map<double,size_t> AgeGroupData::ageMap = fillAgeGroups();

// weight proportions, used by drug code
const double AgeGroupData::wtprop[AgeGroupData::nages] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };


void AgeGroupData::initParameters () {
    for (size_t i=0; i < AgeGroupData::nages; i++) {
	ageSpecificRelativeAvailability[i] = bsa_prop[i] / (1-bsa_prop[i]);
    }
}

const map<double,size_t> AgeGroupData::fillAgeGroups() {
	map<double,size_t> temp;
	for(size_t i=0;i<nages; i++) {
		temp[agemax[i]] = i;
	}
	return temp;
}

size_t AgeGroupData::getAgeGroup (double age) {
    //TODO: use faster search algorithm and compare performance
    /*for (size_t i = 0; i < nages; ++i) {
	if (agemax[i] > age)
	    return i;
    }
    return nages-1;	// final category*/

	map<double, size_t>::const_iterator it;

	it = ageMap.upper_bound(age);

	if(it != ageMap.end()){
		return (*it).second;
	}
	else {
		return nages-1;
	}

}

double AgeGroupData::getAgeSpecificRelativeAvailability (double ageYears) {
    return ageSpecificRelativeAvailability[getAgeGroup (ageYears)];
}

double AgeGroupData::ageToWeight (double ageYears) {

	size_t age_group = getAgeGroup(ageYears);
	size_t lower_group;

	if(age_group>0&&ageYears<60.99)
		lower_group = age_group-1;
	else
		lower_group = age_group;

	double diff_wtprop = wtprop[age_group]-wtprop[lower_group];
	double diff_age_max = agemax[age_group]-agemax[lower_group];
	double diff_age = ageYears - agemax[lower_group];

	double return_value;

	if(diff_age_max>0)
		return_value = 120.0 * (wtprop[lower_group] + ((diff_wtprop/diff_age_max)*diff_age));
	else
		return_value = 120.0 * wtprop[lower_group];

	//double return_value_old = 120 * wtprop[getAgeGroup(ageYears)];

    return return_value;
}
