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


#include "Parameters.h"
#include "util/errors.h"
#include <schema/scenario.h>

#include <boost/format.hpp>

namespace OM {
    using boost::format;
    
// Initialization functions:

Parameters::Parameters( const scnXml::Parameters& parameters ){
    // set parameters
    const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();
    for (scnXml::Parameters::ParameterConstIterator it = paramSeq.begin(); it != paramSeq.end(); ++it) {
        int i = it->getNumber();
        if (i < 0 || i >= MAX)
            continue;   // ignore the parameter; no real point in making this an error
        Parameter parameter = static_cast<Parameter>(i);
        if( !parameterValues.insert( make_pair( parameter, it->getValue() ) ).second )
            throw util::xml_scenario_error( (format("parameter with index %1% described twice") %parameter).str() );
    }
}


double Parameters::operator[]( Parameter parameter )const{
    std::map<Parameter, double>::const_iterator it = parameterValues.find( parameter );
    if( it == parameterValues.end() )
        throw util::xml_scenario_error( (format("parameter %1% required but not described") %parameter).str() );
    return it->second;
}

}
