/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 
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

#include "PkPd/LSTMTreatments.h"
#include "PkPd/LSTMMedicate.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMPkPdModel.h"

#include <boost/foreach.hpp>

namespace OM {
namespace PkPd {

void MedicateData::load( const scnXml::PKPDMedication& med ){
    drug = LSTMDrugType::findDrug( med.getDrug() );
    qty = med.getMg();
    time = med.getHour() / 24.0; // convert from hours to days
    if( med.getDuration().present() ){
        if( !(med.getDuration().get() > 0.0) ){
            throw util::xml_scenario_error( "duration of an IV dose must be some positive amount of time" );
        }
        duration = med.getDuration().get() / 24.0;
    } // else: duration is left as initialised (NaN)
}

void Schedule::load( const scnXml::PKPDSchedule::MedicateSequence& seq ){
    medications.resize( seq.size() );
    size_t i = 0;
    foreach( const scnXml::PKPDMedication& med, seq ){
        medications[i].load( med );
        i += 1;
    }
}

vector<Schedule> schedules;
map<string,size_t> scheduleNames;

void DosageTable::load( const xsd::cxx::tree::sequence<scnXml::PKPDDosageRange>& seq, bool isBodyMass ){
    useMass = isBodyMass;
    double lastMult = 0.0, lastAge = numeric_limits<double>::quiet_NaN();
    foreach( const scnXml::PKPDDosageRange& age, seq ){
        if( lastAge != lastAge ){
            if( age.getLowerbound() != 0.0 )
                throw util::xml_scenario_error( "dosage table must have first lower bound equal 0" );
        }else{
            if( age.getLowerbound() <= lastAge ){
                throw util::xml_scenario_error( "dosage table must list age groups in increasing order" );
            }
            table[age.getLowerbound()] = lastMult;
        }
        lastMult = age.getDose_mult();
        lastAge = age.getLowerbound();
    }
    table[numeric_limits<double>::infinity()] = lastMult;
}

vector<DosageTable> dosages;
map<string,size_t> dosagesNames;

void LSTMTreatments::init(const scnXml::Treatments& data){
    schedules.resize( data.getSchedule().size() );
    size_t i = 0;
    BOOST_FOREACH( const scnXml::PKPDSchedule& schedule, data.getSchedule() ){
        schedules[i].load( schedule.getMedicate() );
        scheduleNames[schedule.getName()] = i;
        i += 1;
    }
    
    dosages.resize( data.getDosages().size() );
    i = 0;
    BOOST_FOREACH( const scnXml::PKPDDosages& elt, data.getDosages() ){
        if( elt.getAge().size() ) dosages[i].load( elt.getAge(), false );
        else dosages[i].load( elt.getBodymass(), true );
        dosagesNames[elt.getName()] = i;
        i += 1;
    }
}

void LSTMTreatments::clear(){
    schedules.clear();
    scheduleNames.clear();
    dosages.clear();
    dosagesNames.clear();
}

size_t LSTMTreatments::findSchedule(const string& name){
    if( scheduleNames.size() == 0 ){
        throw util::xml_scenario_error( "use of PK/PD treatment requires the INCLUDES_PK_PD model option and specification of pharmacology/treatments" );
    }
    map<string,size_t>::const_iterator it = scheduleNames.find( name );
    if( it == scheduleNames.end() ){
        throw util::xml_scenario_error(string("no treatment schedule with this name: ")
            .append(name));
    }
    return it->second;
}

size_t LSTMTreatments::findDosages(const string& name){
    if( dosagesNames.size() == 0 ){
        throw util::xml_scenario_error( "use of PK/PD treatment requires the INCLUDES_PK_PD model option and specification of pharmacology/treatments" );
    }
    map<string,size_t>::const_iterator it = dosagesNames.find( name );
    if( it == dosagesNames.end() ){
        throw util::xml_scenario_error(string("no dosage table with this name: ")
            .append(name));
    }
    return it->second;
}

}
}
