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
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMPkPdModel.h"
#include "util/errors.h"

#include "schema/pharmacology.h"

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
    }
}

struct Schedule {
    void load( const scnXml::PKPDSchedule::MedicateSequence& seq ){
        medications.resize( seq.size() );
        size_t i = 0;
        foreach( const scnXml::PKPDMedication& med, seq ){
            medications[i].load( med );
            i += 1;
        }
    }
    
    vector<MedicateData> medications;
};

vector<Schedule> schedules;
map<string,size_t> scheduleNames;

struct DosageTable {
    void load( const scnXml::PKPDDosages::AgeSequence& seq ){
        double lastMult = 0.0, lastAge = numeric_limits<double>::quiet_NaN();
        foreach( const scnXml::Age& age, seq ){
            if( lastAge != lastAge ){
                if( age.getLb() != 0.0 )
                    throw util::xml_scenario_error( "dosage table must have first lower bound equal 0" );
            }else{
                if( age.getLb() <= lastAge ){
                    throw util::xml_scenario_error( "dosage table must list age groups in increasing order" );
                }
                table[age.getLb()] = lastMult;
            }
            lastMult = age.getDose_mult();
            lastAge = age.getLb();
        }
        table[numeric_limits<double>::infinity()] = lastMult;
    }
    
    double getMultiplier( double age ){
        map<double,double>::const_iterator it = table.upper_bound( age );
        if( it == table.end() )
            throw TRACED_EXCEPTION( "bad age/dosage table", util::Error::PkPd );
        return it->second;
    }
    
    map<double,double> table;
};

vector<DosageTable> dosages;
map<string,size_t> dosagesNames;

void LSTMTreatments::init(const scnXml::Treatments& data)
{
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
        dosages[i].load( elt.getAge() );
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
    map<string,size_t>::const_iterator it = scheduleNames.find( name );
    if( it == scheduleNames.end() ){
        throw util::xml_scenario_error(string("no treatment schedule with this name: ")
            .append(name));
    }
    return it->second;
}

size_t LSTMTreatments::findDosages(const string& name){
    map<string,size_t>::const_iterator it = dosagesNames.find( name );
    if( it == dosagesNames.end() ){
        throw util::xml_scenario_error(string("no dosage table with this name: ")
            .append(name));
    }
    return it->second;
}


void LSTMPkPdModel::prescribe(size_t schedule, size_t dosage, double age){
    double doseMult = dosages[dosage].getMultiplier( age );
    foreach( MedicateData& medicateData, schedules[schedule].medications ){
        medicateQueue.push_back( medicateData.multiplied(doseMult) );
    }
}

}
}
