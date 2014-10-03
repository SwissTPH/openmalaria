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

#include "WHMock.h"
#include "WithinHost/Diagnostic.h"
#include "util/errors.h"

using namespace std;

namespace OM {
namespace UnitTest {

WHMock::WHMock() :
    totalDensity(numeric_limits<double>::quiet_NaN()),
    nTreatments(0)
{}
WHMock::~WHMock() {}

double WHMock::probTransmissionToMosquito( double tbvFactor ) const
{
    throw util::unimplemented_exception( "not needed in unit test" );
}

bool WHMock::summarize(const Host::Human& human){
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::importInfection(){
    throw util::unimplemented_exception( "not needed in unit test" );
}

bool WHMock::optionalPqTreatment(){
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::treatSimple(SimTime timeLiver, SimTime timeBlood){
    nTreatments += 1;
    lastTimeLiver = timeLiver;
    lastTimeBlood = timeBlood;
}

void WHMock::treatPkPd(size_t schedule, size_t dosages, double age){
    nTreatments += 1;
    pkpd.prescribe( schedule, dosages, age );
}

void WHMock::update(int nNewInfs, double ageInYears, double bsvFactor, ofstream& drugMon){
    throw util::unimplemented_exception( "not needed in unit test" );
}

inline double WHMock::getTotalDensity() const{
    return totalDensity;
}

bool WHMock::diagnosticDefault() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::treatment( Host::Human& human, TreatmentId treatId ){
    nTreatments += 1;
}

Pathogenesis::StatePair WHMock::determineMorbidity( double ageYears ){
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::clearImmunity(){
    throw util::unimplemented_exception( "not needed in unit test" );
}

double WHMock::getCumulative_h() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

double WHMock::getCumulative_Y() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

WHInterface::InfectionCount WHMock::countInfections() const{
    throw util::unimplemented_exception( "not needed in unit test" );
}

void WHMock::checkpoint (istream& stream){
    throw util::unimplemented_exception( "not needed in unit test" );
}
void WHMock::checkpoint (ostream& stream){
    throw util::unimplemented_exception( "not needed in unit test" );
}

}
}
